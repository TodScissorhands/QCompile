#include "chemistry/h2_energy.h"
#include "simulator/simulator.h"
#include "simulator/simulation_session.h"
#include <QAbstractItemView>
#include <QApplication>
#include <QButtonGroup>
#include <QDebug>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileInfo>
#include <QFormLayout>
#include <QFont>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMainWindow>
#include <QMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QProgressBar>
#include <QProcess>
#include <QPushButton>
#include <QRadialGradient>
#include <QRadioButton>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QScrollArea>
#include <QSlider>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <QWidget>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <vector>
#include <limits>
#include <iostream>

namespace {

constexpr double kEnergySumTolerance = 5e-12;
constexpr int kProcessOutputLimit = 1 << 20;
constexpr std::array<double, 5> kSweepDistances{0.735, 0.900, 1.100, 1.250, 1.400};

struct EnergyResult {
  double distanceAngstrom = 0.0;
  double distanceBohr = 0.0;
  double parameter = 0.0;
  double electronicEnergy = 0.0;
  double exactElectronicEnergy = 0.0;
  double referenceElectronicEnergy = 0.0;
  double variationalError = 0.0;
  double nuclearRepulsionEnergy = 0.0;
  double totalEnergy = 0.0;
  double exactTotalEnergy = 0.0;
  int evaluations = 0;
  bool converged = false;
};

struct JobToken {
  std::uint64_t id = 0;
  std::atomic<bool> cancelled{false};
  bool sweep = false;
};

struct UiMessage {
  enum class Type { Progress, Point, Complete, Failed, Cancelled } type;
  std::uint64_t id = 0;
  int completed = 0;
  int total = 0;
  bool sweep = false;
  QString text;
  EnergyResult result;
};

QString repositoryRoot() {
  const QString appDirectory = QApplication::applicationDirPath();
  const QStringList candidates{appDirectory, QDir::currentPath(),
                               QDir(appDirectory).filePath("..")};
  for (const QString &candidate : candidates) {
    const QString root = QDir(candidate).absolutePath();
    if (QFileInfo::exists(QDir(root).filePath("Makefile")) &&
        QFileInfo::exists(QDir(root).filePath("benchmarks/h2_geometry/generate_h2.py")) &&
        QFileInfo::exists(QDir(root).filePath("src/apps/benchmark_main.cpp"))) {
      return root;
    }
  }
  throw std::runtime_error("Cannot locate the QCompile repository root; run qcompile_lab from the project tree");
}

QString configuredPath(const char *environmentName, const QString &fallback) {
  const QByteArray value = qgetenv(environmentName);
  return value.isEmpty() ? fallback : QString::fromLocal8Bit(value);
}

struct ProcessResult {
  QByteArray standardOutput;
  QByteArray standardError;
};

ProcessResult runProcess(const QString &stage, const QString &program,
                         const QStringList &arguments,
                         const QString &workingDirectory,
                         const std::shared_ptr<JobToken> &token) {
  QProcess process;
  process.setProgram(program);
  process.setArguments(arguments);
  process.setWorkingDirectory(workingDirectory);
  process.setProcessChannelMode(QProcess::SeparateChannels);
  process.start();
  if (!process.waitForStarted(10000)) {
    throw std::runtime_error(
        QString("%1 process could not start (%2): %3")
            .arg(stage).arg(program).arg(process.errorString()).toStdString());
  }

  ProcessResult result;
  auto capture = [&]() {
    result.standardOutput += process.readAllStandardOutput();
    result.standardError += process.readAllStandardError();
    if (result.standardOutput.size() > kProcessOutputLimit ||
        result.standardError.size() > kProcessOutputLimit) {
      process.kill();
      process.waitForFinished();
      throw std::runtime_error("Child process exceeded the 1 MiB output limit");
    }
  };

  while (process.state() != QProcess::NotRunning) {
    if (token->cancelled.load()) {
      process.kill();
      process.waitForFinished();
      throw std::runtime_error("Calculation cancelled");
    }
    process.waitForFinished(100);
    capture();
  }
  capture();
  if (token->cancelled.load()) {
    throw std::runtime_error("Calculation cancelled");
  }
  if (process.exitStatus() != QProcess::NormalExit || process.exitCode() != 0) {
    const QString details = QString::fromUtf8(result.standardError).trimmed();
    throw std::runtime_error(
        QString("%1 process failed with exit code %2%3")
            .arg(stage)
            .arg(process.exitCode())
            .arg(details.isEmpty() ? QString() : QString(": %1").arg(details))
            .toStdString());
  }
  return result;
}

double requiredJsonNumber(const QJsonObject &object, const QString &name) {
  const QJsonValue value = object.value(name);
  const double number = value.toDouble(std::numeric_limits<double>::quiet_NaN());
  if (!value.isDouble() || !std::isfinite(number)) {
    throw std::runtime_error(QString("JSON result has invalid numeric field: %1")
                                 .arg(name).toStdString());
  }
  return number;
}

void requireJsonString(const QJsonObject &object, const QString &name,
                       const QString &expected) {
  const QJsonValue value = object.value(name);
  if (!value.isString() || value.toString() != expected) {
    throw std::runtime_error(QString("JSON result has incompatible field: %1")
                                 .arg(name).toStdString());
  }
}

EnergyResult parseBenchmarkJson(const QByteArray &output, double requestedDistance) {
  QJsonParseError parseError;
  const QJsonDocument document = QJsonDocument::fromJson(output, &parseError);
  if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
    throw std::runtime_error("Benchmark executable returned malformed JSON");
  }
  const QJsonObject object = document.object();
  requireJsonString(object, QStringLiteral("schema"),
                    QStringLiteral("QCOMPILE-GENERATED-H2-RESULT-1"));
  requireJsonString(object, QStringLiteral("benchmark"), QStringLiteral("generated_h2"));
  requireJsonString(object, QStringLiteral("formula"), QStringLiteral("H2"));
  requireJsonString(object, QStringLiteral("geometry_units"), QStringLiteral("Angstrom"));
  requireJsonString(object, QStringLiteral("basis"), QStringLiteral("STO-3G"));
  requireJsonString(object, QStringLiteral("mapping"), QStringLiteral("parity"));
  requireJsonString(object, QStringLiteral("qubit_ordering"),
                    QStringLiteral("q0 is least-significant state-vector bit"));
  requireJsonString(object, QStringLiteral("energy_unit"), QStringLiteral("Hartree"));
  if (requiredJsonNumber(object, QStringLiteral("charge")) != 0.0 ||
      requiredJsonNumber(object, QStringLiteral("electron_count")) != 2.0 ||
      requiredJsonNumber(object, QStringLiteral("total_spin")) != 0.0 ||
      requiredJsonNumber(object, QStringLiteral("multiplicity")) != 1.0 ||
      requiredJsonNumber(object, QStringLiteral("spatial_orbitals")) != 2.0 ||
      requiredJsonNumber(object, QStringLiteral("qubits")) != 2.0 ||
      !object.value(QStringLiteral("electronic_energy_excludes_nuclear_repulsion")).toBool()) {
    throw std::runtime_error("Benchmark JSON describes an unsupported H2 model");
  }
  const QJsonArray reduction = object.value(QStringLiteral("particle_number_reduction")).toArray();
  if (reduction.size() != 2 || reduction[0].toInt(-1) != 1 || reduction[1].toInt(-1) != 1) {
    throw std::runtime_error("Benchmark JSON has an incompatible parity reduction");
  }

  EnergyResult result;
  result.distanceAngstrom = requiredJsonNumber(object, QStringLiteral("distance_angstrom"));
  result.distanceBohr = requiredJsonNumber(object, QStringLiteral("distance_bohr"));
  if (result.distanceAngstrom <= 0.0 || result.distanceBohr <= 0.0 ||
      std::abs(result.distanceAngstrom - requestedDistance) > 1e-10) {
    throw std::runtime_error("Benchmark geometry does not match the submitted calculation");
  }
  result.parameter = requiredJsonNumber(object, QStringLiteral("optimizer_parameter_radians"));
  result.electronicEnergy = requiredJsonNumber(object, QStringLiteral("electronic_energy_hartree"));
  result.exactElectronicEnergy = requiredJsonNumber(object, QStringLiteral("exact_electronic_energy_hartree"));
  result.referenceElectronicEnergy = requiredJsonNumber(object, QStringLiteral("reference_electronic_energy_hartree"));
  result.variationalError = requiredJsonNumber(object, QStringLiteral("variational_error_hartree"));
  result.nuclearRepulsionEnergy = requiredJsonNumber(object, QStringLiteral("nuclear_repulsion_hartree"));
  result.totalEnergy = requiredJsonNumber(object, QStringLiteral("total_energy_hartree"));
  result.exactTotalEnergy = requiredJsonNumber(object, QStringLiteral("exact_total_energy_hartree"));
  const double evaluationCount = requiredJsonNumber(object, QStringLiteral("evaluation_count"));
  if (evaluationCount < 1.0 || evaluationCount > 10000.0 ||
      std::floor(evaluationCount) != evaluationCount) {
    throw std::runtime_error("Benchmark JSON contains an invalid evaluation count");
  }
  result.evaluations = static_cast<int>(evaluationCount);
  const QJsonValue converged = object.value(QStringLiteral("converged"));
  if (!converged.isBool()) {
    throw std::runtime_error("Benchmark JSON omitted convergence status");
  }
  result.converged = converged.toBool();
  if (!result.converged || result.nuclearRepulsionEnergy <= 0.0 ||
      std::abs(result.totalEnergy - result.electronicEnergy - result.nuclearRepulsionEnergy) > kEnergySumTolerance ||
      std::abs(result.exactTotalEnergy - result.exactElectronicEnergy - result.nuclearRepulsionEnergy) > kEnergySumTolerance ||
      std::abs(result.variationalError - result.electronicEnergy + result.exactElectronicEnergy) > kEnergySumTolerance ||
      std::abs(result.exactElectronicEnergy - result.referenceElectronicEnergy) > 1e-8) {
    throw std::runtime_error("Benchmark result failed energy, reference, or convergence validation");
  }
  return result;
}

struct ProgressCallback {
  std::function<void(const QString &)> publish;
};

EnergyResult calculateAtGeometry(const QString &root, double distance,
                                 const std::shared_ptr<JobToken> &token,
                                 const ProgressCallback &progress) {
  if (!std::isfinite(distance) || distance < 0.4 || distance > 3.0) {
    throw std::invalid_argument("H2 bond length must be in [0.4, 3.0] Angstrom");
  }
  if (token->cancelled.load()) {
    throw std::runtime_error("Calculation cancelled");
  }
  QTemporaryDir temporaryDirectory(QDir::tempPath() + QStringLiteral("/qcompile-h2-lab-XXXXXX"));
  temporaryDirectory.setAutoRemove(true);
  if (!temporaryDirectory.isValid()) {
    throw std::runtime_error("Could not create a private temporary directory for the calculation");
  }
  const QString qhamPath = QDir(temporaryDirectory.path()).filePath("input.qham");
  const QString python = configuredPath(
      "QCOMPILE_LAB_PYTHON",
      QDir(root).filePath("benchmarks/heh_plus/.venv-generator/bin/python"));
  const QString generator = configuredPath(
      "QCOMPILE_LAB_GENERATOR",
      QDir(root).filePath("benchmarks/h2_geometry/generate_h2.py"));
  const QString benchmark = configuredPath(
      "QCOMPILE_LAB_BENCHMARKS",
      QDir(root).filePath("qcompile_benchmarks"));
  if (!QFileInfo(python).isFile() || !QFileInfo(python).isExecutable()) {
    throw std::runtime_error(
        QString("Python interpreter is missing or not executable: %1. Create the pinned generator environment described in benchmarks/heh_plus/GATE_A_STATUS.md from benchmarks/heh_plus/requirements-generator.lock, or set QCOMPILE_LAB_PYTHON to a compatible interpreter.")
            .arg(python).toStdString());
  }
  if (!QFileInfo(generator).isFile()) {
    throw std::runtime_error(
        QString("H2 generator script is missing: %1 (set QCOMPILE_LAB_GENERATOR to override)")
            .arg(generator).toStdString());
  }
  if (!QFileInfo(benchmark).isFile() || !QFileInfo(benchmark).isExecutable()) {
    throw std::runtime_error(
        QString("qcompile_benchmarks executable is missing or not executable: %1 (build it with make qcompile_benchmarks or set QCOMPILE_LAB_BENCHMARKS to override)")
            .arg(benchmark).toStdString());
  }

  progress.publish(QString("Generating STO-3G Hamiltonian at %1 Angstrom…")
                       .arg(distance, 0, 'f', 3));
  runProcess(QStringLiteral("H2 generator"), python,
             {QStringLiteral("-B"), generator, QStringLiteral("--distance-angstrom"),
              QString::number(distance, 'g', 17), QStringLiteral("--output"), qhamPath},
             root, token);
  if (!QFileInfo::exists(qhamPath)) {
    throw std::runtime_error("Generator exited successfully but did not create its Hamiltonian file");
  }
  if (token->cancelled.load()) {
    throw std::runtime_error("Calculation cancelled");
  }
  progress.publish(QString("Running QCompile on generated Hamiltonian at %1 Angstrom…")
                       .arg(distance, 0, 'f', 3));
  const ProcessResult calculation = runProcess(
      QStringLiteral("QCompile benchmark"), benchmark,
      {QStringLiteral("--generated-h2-json"), qhamPath}, root, token);
  return parseBenchmarkJson(calculation.standardOutput, distance);
}

QString formattedEnergy(double value) {
  return QString::number(value, 'f', 12) + QStringLiteral(" Ha");
}

class MoleculeView final : public QWidget {
public:
  explicit MoleculeView(QWidget *parent = nullptr) : QWidget(parent) {
    setMinimumSize(360, 250);
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
  }

  void setDistance(double distance) {
    distance_ = distance;
    update();
  }
  double distance() const { return distance_; }
  double yaw() const { return yaw_; }
  double pitch() const { return pitch_; }
  double zoom() const { return zoom_; }

protected:
  void paintEvent(QPaintEvent *) override {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // Workbench viewport background
    painter.fillRect(rect(), QColor(11, 15, 23));

    // Subtle technical grid pattern
    painter.setPen(QPen(QColor(22, 30, 44), 1, Qt::SolidLine));
    constexpr int kGridStep = 32;
    for (int x = 0; x < width(); x += kGridStep) {
      painter.drawLine(x, 0, x, height());
    }
    for (int y = 0; y < height(); y += kGridStep) {
      painter.drawLine(0, y, width(), y);
    }

    // Viewport border
    painter.setPen(QPen(QColor(30, 41, 59), 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(rect().adjusted(0, 0, -1, -1));

    // Header badge / identification
    painter.setPen(QColor(90, 110, 133));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));
    painter.drawText(QRectF(14, 12, width() - 28, 16), Qt::AlignLeft,
                     QStringLiteral("GEOMETRY SCHEMATIC // H2 MOLECULAR AXIS"));

    const double centerX = width() / 2.0;
    const double centerY = height() * 0.50;
    const double halfDistance = distance_ * 0.5;
    const double cosineYaw = std::cos(yaw_);
    const double sineYaw = std::sin(yaw_);
    const double cosinePitch = std::cos(pitch_);
    const double sinePitch = std::sin(pitch_);
    auto project = [&](double x) {
      const double rotatedX = x * cosineYaw;
      const double rotatedZ = -x * sineYaw;
      const double rotatedY = -rotatedZ * sinePitch;
      const double depth = 3.1 + rotatedZ * cosinePitch;
      const double scale = std::min(width(), height()) * 0.54 * zoom_ / depth;
      return QPointF(centerX + rotatedX * scale, centerY + rotatedY * scale);
    };
    const QPointF first = project(-halfDistance);
    const QPointF second = project(halfDistance);

    // 1. Bond line: clean technical cylindrical axis with centerline
    painter.setPen(QPen(QColor(42, 58, 82), 6, Qt::SolidLine, Qt::RoundCap));
    painter.drawLine(first, second);
    painter.setPen(QPen(QColor(84, 110, 144), 2, Qt::SolidLine, Qt::FlatCap));
    painter.drawLine(first, second);

    // Subtle dashed axis line extending past nuclei
    const QPointF dir = second - first;
    const double len = std::hypot(dir.x(), dir.y());
    if (len > 1e-4) {
      const QPointF unit = dir / len;
      painter.setPen(QPen(QColor(40, 56, 78), 1, Qt::DashLine));
      painter.drawLine(first - unit * 28.0, second + unit * 28.0);
    }

    // 2. Bond dimension readout with calipers
    const QPointF midPoint = (first + second) * 0.5;
    const double perpX = -(second.y() - first.y()) / (len > 1e-4 ? len : 1.0);
    const double perpY = (second.x() - first.x()) / (len > 1e-4 ? len : 1.0);
    const double offset = 26.0;
    const QPointF p1 = first + QPointF(perpX * offset, perpY * offset);
    const QPointF p2 = second + QPointF(perpX * offset, perpY * offset);
    const QPointF pMid = midPoint + QPointF(perpX * offset, perpY * offset);

    painter.setPen(QPen(QColor(56, 189, 248, 160), 1, Qt::DotLine));
    painter.drawLine(first, p1);
    painter.drawLine(second, p2);
    painter.setPen(QPen(QColor(56, 189, 248, 200), 1, Qt::SolidLine));
    painter.drawLine(p1, p2);

    // Caliper ticks
    const QPointF tick(perpX * 3.5, perpY * 3.5);
    painter.drawLine(p1 - tick, p1 + tick);
    painter.drawLine(p2 - tick, p2 + tick);

    // Bond readout badge
    const QString bondStr = QStringLiteral("R = %1 Å").arg(distance_, 0, 'f', 3);
    painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));
    const QFontMetricsF fm(painter.font());
    const double textW = fm.horizontalAdvance(bondStr) + 10.0;
    const double textH = fm.height() + 4.0;
    const QRectF badgeRect(pMid.x() - textW * 0.5, pMid.y() - textH * 0.5, textW, textH);
    painter.setPen(QPen(QColor(56, 189, 248), 1));
    painter.setBrush(QColor(11, 15, 23, 230));
    painter.drawRect(badgeRect);
    painter.setPen(QColor(56, 189, 248));
    painter.drawText(badgeRect, Qt::AlignCenter, bondStr);

    // 3. Atom markers: scientific annular target reticles instead of glossy spheres
    const double atomRadius = std::clamp(std::min(width(), height()) * 0.045, 13.0, 20.0);
    auto drawScientificAtom = [&](const QPointF &center, const QString &label,
                                  const QColor &accentColor) {
      // Outer reticle circle
      painter.setPen(QPen(accentColor.darker(140), 1, Qt::DashLine));
      painter.setBrush(Qt::NoBrush);
      painter.drawEllipse(center, atomRadius + 5.0, atomRadius + 5.0);

      // Core disc: matte technical look
      painter.setPen(QPen(accentColor, 1.5));
      painter.setBrush(QColor(15, 23, 38));
      painter.drawEllipse(center, atomRadius, atomRadius);

      // Reticle tick crosshairs
      painter.setPen(QPen(accentColor, 1));
      painter.drawLine(center.x() - atomRadius - 3, center.y(),
                       center.x() - atomRadius + 3, center.y());
      painter.drawLine(center.x() + atomRadius - 3, center.y(),
                       center.x() + atomRadius + 3, center.y());
      painter.drawLine(center.x(), center.y() - atomRadius - 3,
                       center.x(), center.y() - atomRadius + 3);
      painter.drawLine(center.x(), center.y() + atomRadius - 3,
                       center.x(), center.y() + atomRadius + 3);

      // Label text
      painter.setPen(QColor(241, 245, 249));
      painter.setFont(QFont(QStringLiteral("Monospace"), 9, QFont::Bold));
      const QRectF sphere(center.x() - atomRadius, center.y() - atomRadius,
                          atomRadius * 2.0, atomRadius * 2.0);
      painter.drawText(sphere, Qt::AlignCenter, label);
    };

    drawScientificAtom(first, QStringLiteral("H1"), QColor(56, 189, 248));
    drawScientificAtom(second, QStringLiteral("H2"), QColor(245, 158, 11));

    // Bottom telemetry / scientific disclaimer bar
    painter.setPen(QColor(100, 116, 139));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8));
    painter.drawText(QRectF(14, height() - 32, width() - 28, 14), Qt::AlignLeft,
                     QStringLiteral("ROT: yaw=%1° pitch=%2° · ZOOM: %3x · [LMB drag: rotate / Wheel: zoom]")
                         .arg(yaw_ * 180.0 / 3.1415926535, 5, 'f', 1)
                         .arg(pitch_ * 180.0 / 3.1415926535, 5, 'f', 1)
                         .arg(zoom_, 4, 'f', 2));
    painter.setPen(QColor(71, 85, 105));
    painter.drawText(QRectF(14, height() - 18, width() - 28, 14), Qt::AlignLeft,
                     QStringLiteral("Schematic internuclear vector only; no electron density, orbital, or trajectory model"));
  }

  void mousePressEvent(QMouseEvent *event) override {
    if (event->button() == Qt::LeftButton) {
      dragging_ = true;
      lastPosition_ = event->position();
      event->accept();
      return;
    }
    QWidget::mousePressEvent(event);
  }

  void mouseMoveEvent(QMouseEvent *event) override {
    if (dragging_ && (event->buttons() & Qt::LeftButton)) {
      const QPointF delta = event->position() - lastPosition_;
      yaw_ += delta.x() * 0.01;
      pitch_ = std::clamp(pitch_ + delta.y() * 0.01, -1.3, 1.3);
      lastPosition_ = event->position();
      update();
      event->accept();
      return;
    }
    QWidget::mouseMoveEvent(event);
  }

  void mouseReleaseEvent(QMouseEvent *event) override {
    if (event->button() == Qt::LeftButton) {
      dragging_ = false;
      event->accept();
      return;
    }
    QWidget::mouseReleaseEvent(event);
  }

  void wheelEvent(QWheelEvent *event) override {
    zoom_ = std::clamp(zoom_ * std::pow(1.0015, event->angleDelta().y()), 0.55, 2.2);
    update();
    event->accept();
  }

private:
  double distance_ = 0.735;
  double yaw_ = -0.28;
  double pitch_ = 0.08;
  double zoom_ = 1.0;
  bool dragging_ = false;
  QPointF lastPosition_;
};

class EnergyPlot final : public QWidget {
public:
  explicit EnergyPlot(QWidget *parent = nullptr) : QWidget(parent) {
    setMinimumHeight(235);
  }

  void setResults(const std::vector<EnergyResult> &results) {
    results_ = results;
    update();
  }

protected:
  void paintEvent(QPaintEvent *) override {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // Dark technical background
    painter.fillRect(rect(), QColor(11, 15, 23));

    // Border
    painter.setPen(QPen(QColor(30, 41, 59), 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(rect().adjusted(0, 0, -1, -1));

    // Panel Header
    painter.setPen(QColor(90, 110, 133));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));
    painter.drawText(QRect(12, 10, width() - 24, 16), Qt::AlignLeft,
                     QStringLiteral("BOND-LENGTH ENERGY SCAN (STO-3G) // DISCRETE RUN SAMPLES"));

    if (results_.empty()) {
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 9));
      painter.drawText(rect(), Qt::AlignCenter,
                       QStringLiteral("No calculation samples. Run a geometry or 5-point sweep."));
      return;
    }

    double minX = results_.front().distanceAngstrom;
    double maxX = minX;
    double minY = results_.front().totalEnergy;
    double maxY = minY;
    for (const auto &result : results_) {
      minX = std::min(minX, result.distanceAngstrom);
      maxX = std::max(maxX, result.distanceAngstrom);
      minY = std::min({minY, result.totalEnergy, result.exactTotalEnergy});
      maxY = std::max({maxY, result.totalEnergy, result.exactTotalEnergy});
    }
    if (maxX - minX < 1e-6) { minX -= 0.1; maxX += 0.1; }
    if (maxY - minY < 1e-8) { minY -= 0.01; maxY += 0.01; }
    const double yPad = (maxY - minY) * 0.18;
    minY -= yPad;
    maxY += yPad;

    const QRectF chart(72, 34, std::max(10, width() - 92), std::max(20, height() - 95));

    // Chart background
    painter.fillRect(chart, QColor(15, 23, 42, 160));
    painter.setPen(QPen(QColor(42, 58, 82), 1));
    painter.drawRect(chart);

    // Horizontal grid lines & Y labels
    for (int i = 0; i <= 4; ++i) {
      const double fraction = i / 4.0;
      const double y = chart.bottom() - fraction * chart.height();
      painter.setPen(QPen(QColor(30, 45, 68), 1, Qt::DashLine));
      painter.drawLine(QPointF(chart.left(), y), QPointF(chart.right(), y));
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 8));
      painter.drawText(QRectF(2, y - 7, 66, 14), Qt::AlignRight | Qt::AlignVCenter,
                       QString::number(minY + fraction * (maxY - minY), 'f', 4));
    }

    auto point = [&](double x, double y) {
      return QPointF(chart.left() + (x - minX) / (maxX - minX) * chart.width(),
                     chart.bottom() - (y - minY) / (maxY - minY) * chart.height());
    };

    // Vertical grid lines & X labels
    for (int i = 0; i <= 4; ++i) {
      const double fraction = i / 4.0;
      const double xValue = minX + fraction * (maxX - minX);
      const double xPosition = chart.left() + fraction * chart.width();
      painter.setPen(QPen(QColor(30, 45, 68), 1, Qt::DotLine));
      painter.drawLine(QPointF(xPosition, chart.top()), QPointF(xPosition, chart.bottom()));
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 8));
      painter.drawText(QRectF(xPosition - 25, chart.bottom() + 4, 50, 14),
                       Qt::AlignCenter, QString::number(xValue, 'f', 2));
    }

    // Connect points if more than 1 sample (sorted by distance)
    std::vector<EnergyResult> sortedResults = results_;
    std::sort(sortedResults.begin(), sortedResults.end(),
              [](const EnergyResult &a, const EnergyResult &b) {
                return a.distanceAngstrom < b.distanceAngstrom;
              });

    if (sortedResults.size() > 1) {
      painter.setPen(QPen(QColor(56, 189, 248, 100), 1, Qt::SolidLine));
      for (std::size_t i = 1; i < sortedResults.size(); ++i) {
        painter.drawLine(point(sortedResults[i - 1].distanceAngstrom, sortedResults[i - 1].totalEnergy),
                         point(sortedResults[i].distanceAngstrom, sortedResults[i].totalEnergy));
      }
    }

    const auto minimum = std::min_element(results_.begin(), results_.end(),
        [](const EnergyResult &a, const EnergyResult &b) {
          return a.totalEnergy < b.totalEnergy;
        });

    for (const auto &result : results_) {
      const QPointF calculated = point(result.distanceAngstrom, result.totalEnergy);
      const QPointF reference = point(result.distanceAngstrom, result.exactTotalEnergy);

      // Variational calculated point
      painter.setPen(QPen(QColor(56, 189, 248), 1.5));
      painter.setBrush(QColor(14, 116, 144));
      painter.drawEllipse(calculated, 4.5, 4.5);

      // Reference exact point
      painter.setPen(QPen(QColor(168, 85, 247), 1.5));
      painter.setBrush(Qt::NoBrush);
      QPolygonF triangle;
      triangle << QPointF(reference.x(), reference.y() - 4.5)
               << QPointF(reference.x() - 4.5, reference.y() + 3.5)
               << QPointF(reference.x() + 4.5, reference.y() + 3.5);
      painter.drawPolygon(triangle);
    }

    // Highlight minimum
    if (minimum != results_.end()) {
      painter.setPen(QPen(QColor(245, 158, 11), 1.5));
      painter.setBrush(Qt::NoBrush);
      painter.drawEllipse(point(minimum->distanceAngstrom, minimum->totalEnergy), 8, 8);
    }

    // Chart footer labels & legend
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8));
    painter.drawText(QRectF(chart.left(), chart.bottom() + 18, chart.width(), 16),
                     Qt::AlignCenter, QStringLiteral("BOND LENGTH R (Å)"));

    painter.setPen(QColor(56, 189, 248));
    painter.drawText(QRect(12, height() - 20, 130, 16), Qt::AlignLeft,
                     QStringLiteral("● Variational VQE"));
    painter.setPen(QColor(168, 85, 247));
    painter.drawText(QRect(145, height() - 20, 175, 16), Qt::AlignLeft,
                     QStringLiteral("△ Physical-Sector Exact"));
    if (minimum != results_.end()) {
      painter.setPen(QColor(245, 158, 11));
      painter.drawText(QRect(330, height() - 20, width() - 340, 16), Qt::AlignLeft,
                       QStringLiteral("Sample min: %1 Å (%2 Ha)")
                           .arg(minimum->distanceAngstrom, 0, 'f', 3)
                           .arg(minimum->totalEnergy, 0, 'f', 6));
    }
  }

private:
  std::vector<EnergyResult> results_;
};

class MainWindow final : public QMainWindow {
public:
  explicit MainWindow(bool smokeTest) : root_(repositoryRoot()), smokeTest_(smokeTest) {
    setWindowTitle(QStringLiteral("QCompile — Quantum Chemical Workbench"));
    resize(1440, 920);
    setMinimumSize(1024, 680);
    setStyleSheet(QStringLiteral(
        "QMainWindow, QWidget { background:#0a0e17; color:#cbd5e1; font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; font-size:12px; }"
        "QSplitter::handle { background:#1e293b; }"
        "QSplitter::handle:horizontal { width:2px; }"
        "QSplitter::handle:vertical { height:2px; }"
        "QGroupBox { border:1px solid #1e293b; border-radius:0px; margin-top:8px; padding:10px 8px 8px 8px; font-weight:700; font-size:10px; letter-spacing:0.8px; color:#64748b; }"
        "QGroupBox::title { subcontrol-origin:margin; subcontrol-position:top left; left:8px; padding:0 4px; background:#0a0e17; color:#94a3b8; }"
        "QPushButton { background:#131c2e; border:1px solid #29384e; border-radius:3px; padding:5px 12px; font-weight:600; font-size:11px; color:#e2e8f0; }"
        "QPushButton:hover { background:#1e293b; border-color:#38bdf8; }"
        "QPushButton:pressed { background:#0f172a; }"
        "QPushButton:disabled { background:#0d131f; border-color:#1e293b; color:#475569; }"
        "QPushButton#primary { background:#0284c7; border:1px solid #38bdf8; color:#ffffff; font-weight:700; }"
        "QPushButton#primary:hover { background:#0369a1; border-color:#7dd3fc; }"
        "QPushButton#primary:disabled { background:#1e293b; border-color:#334155; color:#64748b; }"
        "QDoubleSpinBox, QSpinBox { background:#0f172a; border:1px solid #29384e; border-radius:2px; padding:3px 6px; color:#f1f5f9; font-family:monospace; font-size:12px; }"
        "QDoubleSpinBox:focus, QSpinBox:focus { border-color:#38bdf8; }"
        "QTableWidget { background:#0b0f19; border:1px solid #1e293b; border-radius:2px; gridline-color:#1e293b; color:#cbd5e1; font-family:monospace; font-size:11px; selection-background-color:#1e3a5f; selection-color:#f8fafc; }"
        "QHeaderView::section { background:#0f172a; color:#94a3b8; padding:5px 8px; border:0px; border-bottom:1px solid #1e293b; border-right:1px solid #1e293b; font-family:monospace; font-size:10px; font-weight:700; }"
        "QProgressBar { border:1px solid #1e293b; border-radius:2px; text-align:center; background:#0f172a; font-family:monospace; font-size:10px; color:#94a3b8; height:14px; }"
        "QProgressBar::chunk { background:#0284c7; }"
        "QRadioButton { spacing:6px; color:#94a3b8; font-size:11px; }"
        "QRadioButton::indicator { width:12px; height:12px; border:1px solid #334155; border-radius:6px; background:#0f172a; }"
        "QRadioButton::indicator:checked { background:#38bdf8; border-color:#7dd3fc; }"
        "QTabWidget::pane { border:1px solid #1e293b; background:#0b0f19; }"
        "QTabBar::tab { background:#0f172a; border:1px solid #1e293b; border-bottom:0px; padding:6px 14px; color:#64748b; font-weight:600; font-size:10px; letter-spacing:0.5px; }"
        "QTabBar::tab:selected { background:#0b0f19; border-top:2px solid #38bdf8; color:#38bdf8; }"
        "QTabBar::tab:hover:!selected { background:#162032; color:#94a3b8; }"
        "QSlider::groove:horizontal { height:4px; background:#1e293b; border-radius:2px; }"
        "QSlider::sub-page:horizontal { background:#0284c7; border-radius:2px; }"
        "QSlider::handle:horizontal { background:#38bdf8; border:1px solid #7dd3fc; width:12px; margin-top:-4px; margin-bottom:-4px; border-radius:6px; }"
        "QScrollArea { border:0; background:transparent; }"
        "QScrollBar:vertical { background:#0a0e17; width:8px; margin:0px; }"
        "QScrollBar::handle:vertical { background:#1e293b; min-height:20px; border-radius:4px; }"
        "QScrollBar::handle:vertical:hover { background:#334155; }"
        "QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height:0px; }"
        "QScrollBar::horizontal { background:#0a0e17; height:8px; margin:0px; }"
        "QScrollBar::handle:horizontal { background:#1e293b; min-width:20px; border-radius:4px; }"
        "QScrollBar::handle:horizontal:hover { background:#334155; }"
        "QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width:0px; }"));

    buildInterface();
    refreshState();
    refreshCircuit();
    updatePreviewLabels();
    pollTimer_.setInterval(80);
    connect(&pollTimer_, &QTimer::timeout, [this] { drainMessages(); });
    pollTimer_.start();
    playbackTimer_.setInterval(500);
    connect(&playbackTimer_, &QTimer::timeout, [this] { stepForward(); });
    if (smokeTest_) {
      QTimer::singleShot(500, [this] { beginSmokeTest(); });
    }
  }

  ~MainWindow() override {
    if (token_) token_->cancelled.store(true);
    if (worker_.joinable()) worker_.join();
  }

private:
  void buildInterface() {
    auto *central = new QWidget;
    auto *outer = new QVBoxLayout(central);
    outer->setContentsMargins(10, 8, 10, 8);
    outer->setSpacing(6);

    // 1. Compact professional header toolbar
    auto *header = new QHBoxLayout;
    header->setSpacing(12);

    auto *brandLabel = new QLabel(QStringLiteral("QCOMPILE"));
    brandLabel->setStyleSheet(QStringLiteral("font-family:monospace;font-size:14px;font-weight:900;letter-spacing:1.5px;color:#38bdf8"));

    auto *sep1 = new QLabel(QStringLiteral("/"));
    sep1->setStyleSheet(QStringLiteral("color:#334155;font-weight:bold"));

    auto *modLabel = new QLabel(QStringLiteral("H2 MOLECULAR WORKBENCH"));
    modLabel->setStyleSheet(QStringLiteral("font-family:monospace;font-size:12px;font-weight:700;letter-spacing:1px;color:#f1f5f9"));

    auto *sep2 = new QLabel(QStringLiteral("·"));
    sep2->setStyleSheet(QStringLiteral("color:#334155"));

    auto *subLabel = new QLabel(QStringLiteral("STO-3G // Parity Mapping // 2-Qubit Statevector"));
    subLabel->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b"));

    header->addWidget(brandLabel);
    header->addWidget(sep1);
    header->addWidget(modLabel);
    header->addWidget(sep2);
    header->addWidget(subLabel);
    header->addStretch(1);

    appStatus_ = new QLabel(QStringLiteral("READY · STO-3G"));
    appStatus_->setStyleSheet(QStringLiteral(
        "font-family:monospace;font-size:10px;font-weight:700;color:#38bdf8;padding:3px 8px;background:#0f172a;border:1px solid #1e293b;border-radius:2px"));
    header->addWidget(appStatus_);

    outer->addLayout(header);

    // Thin separator line
    auto *line = new QFrame;
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Plain);
    line->setStyleSheet(QStringLiteral("color:#1e293b;max-height:1px"));
    outer->addWidget(line);

    // 2. Main layout: Primary QSplitter (Left: Scientific Viewport & Circuit; Right: Controls & Results)
    auto *columns = new QSplitter(Qt::Horizontal);
    columns->setChildrenCollapsible(false);

    // LEFT PANE: Scientific Visualization & Circuit Inspector
    auto *leftSplitter = new QSplitter(Qt::Vertical);
    leftSplitter->setChildrenCollapsible(false);

    // Upper Left: 3D Molecule Schematic Viewport
    moleculeView_ = new MoleculeView;
    moleculeView_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    leftSplitter->addWidget(moleculeView_);

    // Lower Left: Circuit & State Inspection Pane
    auto *circuitBox = new QWidget;
    auto *circuitLayout = new QVBoxLayout(circuitBox);
    circuitLayout->setContentsMargins(0, 4, 0, 0);
    circuitLayout->setSpacing(6);

    // Composer-inspired compact playback control bar
    auto *playbackBar = new QHBoxLayout;
    playbackBar->setSpacing(6);

    auto *playbackTitle = new QLabel(QStringLiteral("EXECUTION STEP:"));
    playbackTitle->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;font-weight:700;color:#64748b"));
    playbackBar->addWidget(playbackTitle);

    resetButton_ = new QPushButton(QStringLiteral("Reset"));
    stepButton_ = new QPushButton(QStringLiteral("Step ▶|"));
    playButton_ = new QPushButton(QStringLiteral("Play ▶"));
    pauseButton_ = new QPushButton(QStringLiteral("Pause ⏸"));
    resetButton_->setFixedHeight(24);
    stepButton_->setFixedHeight(24);
    playButton_->setFixedHeight(24);
    pauseButton_->setFixedHeight(24);

    playbackBar->addWidget(resetButton_);
    playbackBar->addWidget(stepButton_);
    playbackBar->addWidget(playButton_);
    playbackBar->addWidget(pauseButton_);

    stepSlider_ = new QSlider(Qt::Horizontal);
    stepSlider_->setRange(0, 8);
    stepSlider_->setFixedHeight(18);
    playbackBar->addWidget(stepSlider_, 1);

    circuitLayout->addLayout(playbackBar);

    // Composer tabs: Gate Sequence vs State Inspector
    auto *tabs = new QTabWidget;

    // Gate sequence tab
    auto *gatePage = new QWidget;
    auto *gateLayout = new QVBoxLayout(gatePage);
    gateLayout->setContentsMargins(4, 4, 4, 4);
    gateLayout->setSpacing(4);

    gateTable_ = new QTableWidget(0, 4);
    gateTable_->setHorizontalHeaderLabels({"#", "Instruction", "Qubits", "Parameter"});
    gateTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    gateTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    gateTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    gateTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
    gateTable_->verticalHeader()->hide();
    gateTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    gateTable_->setSelectionMode(QAbstractItemView::NoSelection);
    gateLayout->addWidget(gateTable_);
    tabs->addTab(gatePage, QStringLiteral("Circuit Instructions (Trial Sequence)"));

    // State inspector tab
    auto *statePage = new QWidget;
    auto *stateLayout = new QVBoxLayout(statePage);
    stateLayout->setContentsMargins(4, 4, 4, 4);
    stateLayout->setSpacing(4);

    stateStageLabel_ = new QLabel;
    stateStageLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#38bdf8;padding:2px 0;font-weight:600"));
    stateLayout->addWidget(stateStageLabel_);

    stateTable_ = new QTableWidget(4, 3);
    stateTable_->setHorizontalHeaderLabels({"Basis (q1q0; q0 LSB)", "Amplitude Re + i Im", "Probability |c|²"});
    stateTable_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    stateTable_->verticalHeader()->hide();
    stateTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    stateTable_->setSelectionMode(QAbstractItemView::NoSelection);
    stateLayout->addWidget(stateTable_);

    stateNormLabel_ = new QLabel;
    stateNormLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#94a3b8;padding:2px 0"));
    stateLayout->addWidget(stateNormLabel_);

    tabs->addTab(statePage, QStringLiteral("State Inspector (Statevector)"));
    circuitLayout->addWidget(tabs);

    connect(resetButton_, &QPushButton::clicked, [this] { playbackTimer_.stop(); stepSlider_->setValue(0); });
    connect(stepButton_, &QPushButton::clicked, [this] { playbackTimer_.stop(); stepSlider_->setValue(std::min(8, stepSlider_->value() + 1)); });
    connect(playButton_, &QPushButton::clicked, [this] {
      if (stepSlider_->value() >= 8) stepSlider_->setValue(0);
      playbackTimer_.start();
    });
    connect(pauseButton_, &QPushButton::clicked, [this] { playbackTimer_.stop(); });
    connect(stepSlider_, &QSlider::valueChanged, [this] { refreshState(); refreshCircuit(); });

    leftSplitter->addWidget(circuitBox);
    leftSplitter->setStretchFactor(0, 11);
    leftSplitter->setStretchFactor(1, 9);

    columns->addWidget(leftSplitter);

    // RIGHT PANE: Scientific Controls, Telemetry & Results
    // Right pane tabs: Focus calculation workflow into clean workspaces
    auto *rightTabs = new QTabWidget;

    // Tab 1: Experiment & Results (Geometry control + Latest calculation results)
    auto *calcPage = new QWidget;
    auto *calcLayout = new QVBoxLayout(calcPage);
    calcLayout->setContentsMargins(8, 8, 8, 8);
    calcLayout->setSpacing(8);

    // Geometry & Action Controls Group
    auto *geometry = new QGroupBox(QStringLiteral("Geometry & Calculation Dispatch"));
    auto *geometryLayout = new QVBoxLayout(geometry);
    geometryLayout->setSpacing(6);

    auto *paramRow = new QHBoxLayout;
    auto *distTitle = new QLabel(QStringLiteral("Target Bond Length:"));
    distTitle->setStyleSheet(QStringLiteral("font-family:monospace;font-size:11px;color:#94a3b8"));
    distanceSpin_ = new QDoubleSpinBox;
    distanceSpin_->setRange(0.400, 3.000);
    distanceSpin_->setDecimals(3);
    distanceSpin_->setSingleStep(0.025);
    distanceSpin_->setValue(0.735);
    distanceSpin_->setSuffix(QStringLiteral(" Å"));
    distanceSpin_->setMinimumWidth(110);
    paramRow->addWidget(distTitle);
    paramRow->addWidget(distanceSpin_);
    paramRow->addStretch(1);
    geometryLayout->addLayout(paramRow);

    previewLabel_ = new QLabel;
    previewLabel_->setWordWrap(true);
    previewLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#38bdf8"));
    geometryLayout->addWidget(previewLabel_);

    connect(distanceSpin_, qOverload<double>(&QDoubleSpinBox::valueChanged), [this](double distance) {
      moleculeView_->setDistance(distance);
      updatePreviewLabels();
    });

    auto *btnRow = new QHBoxLayout;
    btnRow->setSpacing(6);
    runButton_ = new QPushButton(QStringLiteral("Calculate Geometry"));
    runButton_->setObjectName(QStringLiteral("primary"));
    runButton_->setFixedHeight(28);

    sweepButton_ = new QPushButton(QStringLiteral("Run 5-Point Sweep"));
    sweepButton_->setFixedHeight(28);

    cancelButton_ = new QPushButton(QStringLiteral("Cancel"));
    cancelButton_->setFixedHeight(28);
    cancelButton_->setEnabled(false);

    btnRow->addWidget(runButton_, 2);
    btnRow->addWidget(sweepButton_, 2);
    btnRow->addWidget(cancelButton_, 1);
    geometryLayout->addLayout(btnRow);

    // Progress bar and status
    auto *progRow = new QHBoxLayout;
    progress_ = new QProgressBar;
    progress_->setRange(0, 100);
    progress_->setTextVisible(true);
    progRow->addWidget(progress_);
    geometryLayout->addLayout(progRow);

    progressLabel_ = new QLabel(QStringLiteral("Ready. Isolated temporary Hamiltonian generated per run."));
    progressLabel_->setWordWrap(true);
    progressLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b"));
    geometryLayout->addWidget(progressLabel_);

    connect(runButton_, &QPushButton::clicked, [this] { startJob(false); });
    connect(sweepButton_, &QPushButton::clicked, [this] { startJob(true); });
    connect(cancelButton_, &QPushButton::clicked, [this] { if (token_) token_->cancelled.store(true); });

    calcLayout->addWidget(geometry);

    // Completed Calculation Metrics Table / Form
    auto *energyBox = new QGroupBox(QStringLiteral("Latest Calculation Results // Energy Breakdown"));
    auto *energyLayout = new QFormLayout(energyBox);
    energyLayout->setLabelAlignment(Qt::AlignLeft);
    energyLayout->setFormAlignment(Qt::AlignLeft | Qt::AlignTop);
    energyLayout->setHorizontalSpacing(16);
    energyLayout->setVerticalSpacing(5);

    energyGeometry_ = addValueRow(energyLayout, QStringLiteral("Calculated Geometry"));
    electronicValue_ = addValueRow(energyLayout, QStringLiteral("Electronic Energy"));
    nuclearValue_ = addValueRow(energyLayout, QStringLiteral("Nuclear Repulsion"));
    totalValue_ = addValueRow(energyLayout, QStringLiteral("Total Born–Oppenheimer"));
    exactValue_ = addValueRow(energyLayout, QStringLiteral("Physical-Sector Exact"));
    referenceValue_ = addValueRow(energyLayout, QStringLiteral("Generator FCI Reference"));
    errorValue_ = addValueRow(energyLayout, QStringLiteral("Variational Error (ΔE)"));
    optimizerValue_ = addValueRow(energyLayout, QStringLiteral("Optimizer / State"));

    calcLayout->addWidget(energyBox);

    // Guided Prediction Group
    auto *question = new QGroupBox(QStringLiteral("Guided Prediction // Energy Order"));
    auto *questionLayout = new QVBoxLayout(question);
    questionLayout->setSpacing(4);

    auto *questionText = new QLabel(QStringLiteral("Hypothesis: compare total Born–Oppenheimer energy across the test interval:"));
    questionText->setWordWrap(true);
    questionText->setStyleSheet(QStringLiteral("font-size:11px;color:#94a3b8"));
    questionLayout->addWidget(questionText);

    auto *radioRow = new QHBoxLayout;
    predictionGroup_ = new QButtonGroup(this);
    predictionLower_ = new QRadioButton(QStringLiteral("E(0.735) < E(1.400)"));
    predictionLonger_ = new QRadioButton(QStringLiteral("E(1.400) < E(0.735)"));
    predictionEqual_ = new QRadioButton(QStringLiteral("E(0.735) ≈ E(1.400)"));
    predictionGroup_->addButton(predictionLower_, 0);
    predictionGroup_->addButton(predictionLonger_, 1);
    predictionGroup_->addButton(predictionEqual_, 2);
    predictionLower_->setChecked(true);
    radioRow->addWidget(predictionLower_);
    radioRow->addWidget(predictionLonger_);
    radioRow->addWidget(predictionEqual_);
    questionLayout->addLayout(radioRow);

    predictionFeedback_ = new QLabel(QStringLiteral("Prediction recorded; validated upon running individual geometries or sweep."));
    predictionFeedback_->setWordWrap(true);
    predictionFeedback_->setStyleSheet(QStringLiteral("font-size:11px;color:#94a3b8;padding:6px 8px;background:#090d16;border:1px solid #1e293b;border-left:2px solid #38bdf8;min-height:48px"));
    questionLayout->addWidget(predictionFeedback_);

    calcLayout->addWidget(question);
    calcLayout->addStretch(1);
    rightTabs_ = rightTabs;
    rightTabs_->addTab(calcPage, QStringLiteral("Calculation & Telemetry"));

    // Tab 2: Energy Scan Plot & Analysis
    auto *plotContainer = new QWidget;
    auto *plotContainerLayout = new QVBoxLayout(plotContainer);
    plotContainerLayout->setContentsMargins(8, 8, 8, 8);
    plotContainerLayout->setSpacing(6);

    energyPlot_ = new EnergyPlot;
    energyPlot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    plotContainerLayout->addWidget(energyPlot_, 1);

    explanation_ = new QLabel(QStringLiteral("Variational energy is compared with the physical-sector exact ground reference of the generated STO-3G Hamiltonian. No quantum advantage is claimed."));
    explanation_->setWordWrap(true);
    explanation_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b;padding:6px 8px;border:1px solid #1e293b;background:#090d16"));
    plotContainerLayout->addWidget(explanation_);

    rightTabs->addTab(plotContainer, QStringLiteral("Energy Scan (5-Point Sweep)"));

    columns->addWidget(rightTabs);

    columns->setStretchFactor(0, 11);
    columns->setStretchFactor(1, 10);

    outer->addWidget(columns, 1);
    setCentralWidget(central);
  }

  QLabel *addValueRow(QFormLayout *layout, const QString &name) {
    auto *nameLabel = new QLabel(name);
    nameLabel->setStyleSheet(QStringLiteral("font-size:11px;color:#94a3b8"));
    auto *value = new QLabel(QStringLiteral("—"));
    value->setTextInteractionFlags(Qt::TextSelectableByMouse);
    value->setStyleSheet(QStringLiteral("font-family:monospace;font-size:11px;color:#f1f5f9;font-weight:600"));
    layout->addRow(nameLabel, value);
    return value;
  }

  void refreshCircuit() {
    const double parameter = latestResult_ ? latestResult_->parameter : 0.0;
    const auto instructions = h2TrialInstructions(parameter);
    gateTable_->setRowCount(static_cast<int>(instructions.size()));
    for (std::size_t index = 0; index < instructions.size(); ++index) {
      const auto &instruction = instructions[index];
      QString qubits;
      for (std::size_t i = 0; i < instruction.qubitIndices.size(); ++i) {
        if (i != 0) qubits += QStringLiteral(", ");
        qubits += QStringLiteral("q%1").arg(instruction.qubitIndices[i]);
      }
      QString parameterText = QStringLiteral("—");
      if (instruction.gateName == "RZ") {
        parameterText = QString::number(instruction.parameterValue, 'f', 9) + QStringLiteral(" rad");
        if (index == 4 && latestResult_) parameterText = QStringLiteral("θ = %1 rad").arg(latestResult_->parameter, 0, 'f', 9);
      }
      const std::array<QString, 4> cells{
          QString::number(index + 1), QString::fromStdString(instruction.gateName), qubits, parameterText};
      for (int column = 0; column < 4; ++column) {
        auto *item = new QTableWidgetItem(cells[static_cast<std::size_t>(column)]);
        item->setTextAlignment(column == 1 ? Qt::AlignLeft | Qt::AlignVCenter : Qt::AlignCenter);
        if (static_cast<int>(index) == stepSlider_->value() && stepSlider_->value() < static_cast<int>(instructions.size())) {
          item->setBackground(QColor(14, 116, 144));
          item->setForeground(QColor(240, 253, 250));
        }
        gateTable_->setItem(static_cast<int>(index), column, item);
      }
    }
  }

  void refreshState() {
    const double parameter = latestResult_ ? latestResult_->parameter : 0.0;
    const auto instructions = h2TrialInstructions(parameter);
    const int completed = stepSlider_->value();
    SimulationSession session(instructions, 2);
    session.seek(static_cast<std::size_t>(completed));
    const auto &state = session.currentState();
    const double norm = session.stateNorm();
    stateStageLabel_->setText(
        completed == 0 ? QStringLiteral("Selected stage: initial |q1q0⟩ = |00⟩ basis state before gate 1")
                       : QStringLiteral("Selected stage: after gate %1 of %2")
                             .arg(completed).arg(instructions.size()));
    for (int index = 0; index < 4; ++index) {
      const auto &amplitude = state[static_cast<std::size_t>(index)];
      const QString label = QStringLiteral("|%1⟩ (q1q0; index %2)")
                                .arg(QString::number(index, 2).rightJustified(2, QLatin1Char('0')))
                                .arg(index);
      stateTable_->setItem(index, 0, new QTableWidgetItem(label));
      stateTable_->setItem(index, 1, new QTableWidgetItem(
          QStringLiteral("%1 %2 %3 i").arg(amplitude.real(), 0, 'f', 9)
              .arg(amplitude.imag() < 0 ? QStringLiteral("−") : QStringLiteral("+"))
              .arg(std::abs(amplitude.imag()), 0, 'f', 9)));
      stateTable_->setItem(index, 2, new QTableWidgetItem(
          QString::number(std::norm(amplitude), 'f', 9)));
    }
    stateNormLabel_->setText(QStringLiteral("Probability = |amplitude|²; norm = %1 · %2")
                                 .arg(norm, 0, 'f', 12)
                                 .arg(std::abs(norm - 1.0) < 1e-12 ? QStringLiteral("normalized") : QStringLiteral("invalid")));
    refreshCircuit();
  }

  void stepForward() {
    if (stepSlider_->value() >= 8) {
      playbackTimer_.stop();
      return;
    }
    stepSlider_->setValue(stepSlider_->value() + 1);
    if (stepSlider_->value() >= 8) playbackTimer_.stop();
  }

  void updatePreviewLabels() {
    if (!moleculeView_) return;
    const QString last = latestResult_
        ? QStringLiteral("Preview %1 Å; latest completed calculation %2 Å")
              .arg(distanceSpin_->value(), 0, 'f', 3)
              .arg(latestResult_->distanceAngstrom, 0, 'f', 3)
        : QStringLiteral("Preview %1 Å; no completed calculation yet")
              .arg(distanceSpin_->value(), 0, 'f', 3);
    previewLabel_->setText(last);
  }

  void post(UiMessage message) {
    std::lock_guard<std::mutex> lock(messagesMutex_);
    messages_.push_back(std::move(message));
  }

  void startJob(bool sweep) {
    if (jobRunning_) return;
    if (worker_.joinable()) worker_.join();
    playbackTimer_.stop();
    auto token = std::make_shared<JobToken>();
    token->id = ++requestId_;
    token->sweep = sweep;
    token_ = token;
    jobRunning_ = true;
    runButton_->setEnabled(false);
    if (sweep) currentSweepResults_.clear();
    sweepButton_->setEnabled(false);
    cancelButton_->setEnabled(true);
    progress_->setValue(0);
    progressLabel_->setText(sweep ? QStringLiteral("Starting fresh five-geometry generator/solver sweep…")
                                  : QStringLiteral("Starting fresh geometry-specific calculation…"));
    appStatus_->setText(QStringLiteral("Calculation running · request %1").arg(token->id));
    const QString root = root_;
    const double selectedDistance = distanceSpin_->value();
    worker_ = std::thread([this, token, sweep, selectedDistance, root] {
      const int total = sweep ? static_cast<int>(kSweepDistances.size()) : 1;
      try {
        for (int index = 0; index < total; ++index) {
          const double distance = sweep ? kSweepDistances[static_cast<std::size_t>(index)] : selectedDistance;
          post({UiMessage::Type::Progress, token->id, index, total, sweep,
                QStringLiteral("Preparing sample %1/%2 at %3 Å")
                    .arg(index + 1).arg(total).arg(distance, 0, 'f', 3), {}});
          const ProgressCallback callback{[this, token, index, total, sweep](const QString &text) {
            post({UiMessage::Type::Progress, token->id, index, total, sweep, text, {}});
          }};
          EnergyResult result = calculateAtGeometry(root, distance, token, callback);
          if (token->cancelled.load()) throw std::runtime_error("Calculation cancelled");
          post({UiMessage::Type::Point, token->id, index + 1, total, sweep,
                QStringLiteral("Completed actual calculation at %1 Å")
                    .arg(result.distanceAngstrom, 0, 'f', 3), result});
        }
        post({UiMessage::Type::Complete, token->id, total, total, sweep,
              QStringLiteral("All requested calculations completed."), {}});
      } catch (const std::exception &error) {
        const bool cancelled = token->cancelled.load() ||
            QString::fromUtf8(error.what()).contains(QStringLiteral("cancelled"), Qt::CaseInsensitive);
        post({cancelled ? UiMessage::Type::Cancelled : UiMessage::Type::Failed,
              token->id, 0, 0, sweep, QString::fromUtf8(error.what()), {}});
      }
    });
  }

  void drainMessages() {
    std::deque<UiMessage> messages;
    {
      std::lock_guard<std::mutex> lock(messagesMutex_);
      messages.swap(messages_);
    }
    for (const UiMessage &message : messages) {
      if (message.id != requestId_) continue;
      switch (message.type) {
      case UiMessage::Type::Progress:
        progressLabel_->setText(message.text);
        if (message.total > 0) progress_->setValue(message.completed * 100 / message.total);
        break;
      case UiMessage::Type::Point:
        plotResults_.push_back(message.result);
        if (message.sweep) currentSweepResults_.push_back(message.result);
        energyPlot_->setResults(plotResults_);
        latestResult_ = message.result;
        std::cout
            << QStringLiteral("QCOMPILE_LAB_RESULT R=%1 A Eel=%2 Ha Enuc=%3 Ha Etot=%4 Ha Eexact=%5 Ha dE=%6 Ha evals=%7 converged=%8")
                   .arg(message.result.distanceAngstrom, 0, 'f', 6)
                   .arg(message.result.electronicEnergy, 0, 'f', 12)
                   .arg(message.result.nuclearRepulsionEnergy, 0, 'f', 12)
                   .arg(message.result.totalEnergy, 0, 'f', 12)
                   .arg(message.result.exactElectronicEnergy, 0, 'f', 12)
                   .arg(message.result.variationalError, 0, 'f', 12)
                   .arg(message.result.evaluations)
                   .arg(message.result.converged ? QStringLiteral("yes") : QStringLiteral("no"))
                   .toStdString() << std::endl;
        updateResultFields(message.result);
        refreshCircuit();
        refreshState();
        updatePreviewLabels();
        progress_->setValue(message.completed * 100 / std::max(1, message.total));
        progressLabel_->setText(message.text);
        appStatus_->setText(QStringLiteral("Completed %1 · run %2")
                                .arg(message.sweep ? QStringLiteral("sweep sample") : QStringLiteral("calculation"))
                                .arg(message.id));
        break;
      case UiMessage::Type::Complete:
        jobRunning_ = false;
        runButton_->setEnabled(true);
        sweepButton_->setEnabled(true);
        cancelButton_->setEnabled(false);
        progress_->setValue(100);
        progressLabel_->setText(message.text);
        if (worker_.joinable()) worker_.join();
        if (message.sweep) {
          updatePredictionFeedback();
          if (rightTabs_) rightTabs_->setCurrentIndex(1);
        }
        if (smokeTest_) continueSmokeTest();
        break;
      case UiMessage::Type::Failed:
      case UiMessage::Type::Cancelled:
        jobRunning_ = false;
        runButton_->setEnabled(true);
        sweepButton_->setEnabled(true);
        cancelButton_->setEnabled(false);
        progressLabel_->setText((message.type == UiMessage::Type::Failed ? QStringLiteral("Failed: ") : QStringLiteral("Cancelled: ")) + message.text);
        appStatus_->setText(message.type == UiMessage::Type::Failed ? QStringLiteral("Run failed; previous result remains labelled by its own geometry") : QStringLiteral("Run cancelled"));
        if (worker_.joinable()) worker_.join();
        if (smokeTest_) finishSmokeTest(false, message.text);
        break;
      }
    }
  }

  void updateResultFields(const EnergyResult &result) {
    energyGeometry_->setText(QStringLiteral("%1 Å (%2 bohr)")
                                 .arg(result.distanceAngstrom, 0, 'f', 6)
                                 .arg(result.distanceBohr, 0, 'f', 9));
    electronicValue_->setText(formattedEnergy(result.electronicEnergy));
    nuclearValue_->setText(formattedEnergy(result.nuclearRepulsionEnergy));
    totalValue_->setText(formattedEnergy(result.totalEnergy));
    exactValue_->setText(formattedEnergy(result.exactElectronicEnergy));
    referenceValue_->setText(formattedEnergy(result.referenceElectronicEnergy));
    errorValue_->setText(formattedEnergy(result.variationalError));
    optimizerValue_->setText(QStringLiteral("%1 evaluations · %2 · θ=%3 rad")
                                 .arg(result.evaluations)
                                 .arg(result.converged ? QStringLiteral("converged") : QStringLiteral("not converged"))
                                 .arg(result.parameter, 0, 'f', 9));
  }

  void updatePredictionFeedback() {
    auto findAt = [&](double distance) -> std::optional<EnergyResult> {
      for (auto it = currentSweepResults_.rbegin(); it != currentSweepResults_.rend(); ++it) {
        if (std::abs(it->distanceAngstrom - distance) < 1e-8) return *it;
      }
      return std::nullopt;
    };
    const auto shortBond = findAt(0.735);
    const auto longBond = findAt(1.400);
    if (!shortBond || !longBond) return;
    const double difference = longBond->totalEnergy - shortBond->totalEnergy;
    const int prediction = predictionGroup_->checkedId();
    const bool correct = prediction == (difference > 1e-6 ? 0 : (difference < -1e-6 ? 1 : 2));
    const auto minimum = std::min_element(currentSweepResults_.begin(), currentSweepResults_.end(),
        [](const EnergyResult &a, const EnergyResult &b) { return a.totalEnergy < b.totalEnergy; });
    predictionFeedback_->setText(
        QStringLiteral("Measured: E(1.400 Å) − E(0.735 Å) = %1 Ha. Your prediction was %2. Sample minimum: %3 Å (%4 Ha). These are five STO-3G calculations; they do not establish accuracy of the one-parameter ansatz throughout dissociation.")
            .arg(difference, 0, 'f', 9)
            .arg(correct ? QStringLiteral("consistent with the sampled endpoints") : QStringLiteral("not consistent with the sampled endpoints"))
            .arg(minimum->distanceAngstrom, 0, 'f', 3)
            .arg(minimum->totalEnergy, 0, 'f', 9));
    explanation_->setText(predictionFeedback_->text() +
        QStringLiteral(" The exact reference is the loaded Hamiltonian’s physical-sector minimum, separately reported from the variational result. No quantum-advantage claim is made."));
  }

  void beginSmokeTest() {
    if (!smokeTest_) return;
    // Exercise the actual camera event handlers and the circuit playback controls.
    const double initialYaw = moleculeView_->yaw();
    const QPoint pressPosition(40, 40);
    QMouseEvent press(QEvent::MouseButtonPress, QPointF(pressPosition),
                      QPointF(moleculeView_->mapToGlobal(pressPosition)),
                      Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(moleculeView_, &press);
    const QPoint movePosition(85, 60);
    QMouseEvent move(QEvent::MouseMove, QPointF(movePosition),
                     QPointF(moleculeView_->mapToGlobal(movePosition)),
                     Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(moleculeView_, &move);
    QMouseEvent release(QEvent::MouseButtonRelease, QPointF(movePosition),
                        QPointF(moleculeView_->mapToGlobal(movePosition)),
                        Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(moleculeView_, &release);
    if (std::abs(moleculeView_->yaw() - initialYaw) < 0.1) {
      finishSmokeTest(false, QStringLiteral("camera drag did not rotate the 3D view"));
      return;
    }
    const double initialZoom = moleculeView_->zoom();
    QWheelEvent wheel(QPointF(100, 100),
                      QPointF(moleculeView_->mapToGlobal(QPoint(100, 100))),
                      QPoint(), QPoint(0, 120), Qt::NoButton, Qt::NoModifier,
                      Qt::NoScrollPhase, false);
    QApplication::sendEvent(moleculeView_, &wheel);
    if (std::abs(moleculeView_->zoom() - initialZoom) < 1e-6) {
      finishSmokeTest(false, QStringLiteral("wheel control did not zoom the molecule view"));
      return;
    }
    stepButton_->click();
    if (stepSlider_->value() != 1) {
      finishSmokeTest(false, QStringLiteral("single-step control did not advance the circuit"));
      return;
    }
    if (!stateTable_->item(1, 2) ||
        std::abs(stateTable_->item(1, 2)->text().toDouble() - 1.0) > 1e-12 ||
        !stateNormLabel_->text().contains(QStringLiteral("normalized"))) {
      finishSmokeTest(false, QStringLiteral("one-gate state inspector did not show the simulator's normalized |01> state"));
      return;
    }
    resetButton_->click();
    playButton_->click();
    QTimer::singleShot(550, [this] {
      pauseButton_->click();
      if (stepSlider_->value() != 1) {
        finishSmokeTest(false, QStringLiteral("play/pause controls did not advance one gate"));
        return;
      }
      distanceSpin_->setValue(0.735);
      if (std::abs(moleculeView_->distance() - 0.735) > 1e-12) {
        finishSmokeTest(false, QStringLiteral("bond-length control did not update the molecule preview"));
        return;
      }
      runButton_->click();
    });
  }

  void continueSmokeTest() {
    ++smokeStage_;
    if (smokeStage_ == 1) {
      smokeFirstResult_ = latestResult_;
      if (!smokeFirstResult_ || std::abs(smokeFirstResult_->distanceAngstrom - 0.735) > 1e-10) {
        finishSmokeTest(false, QStringLiteral("first fresh geometry result was not 0.735 Angstrom"));
        return;
      }
      distanceSpin_->setValue(1.400);
      runButton_->click();
    } else if (smokeStage_ == 2) {
      if (!latestResult_ || !smokeFirstResult_ ||
          std::abs(latestResult_->distanceAngstrom - 1.400) > 1e-10 ||
          std::abs(latestResult_->totalEnergy - smokeFirstResult_->totalEnergy) < 0.1) {
        finishSmokeTest(false, QStringLiteral("distinct 1.400 Angstrom result was not validated"));
        return;
      }
      sweepButton_->click();
    } else {
      if (currentSweepResults_.size() != kSweepDistances.size()) {
        finishSmokeTest(false, QStringLiteral("five-point sweep did not produce five actual results"));
        return;
      }
      finishSmokeTest(true, QStringLiteral("Two distinct geometries and five fresh sweep calculations completed"));
    }
  }

  void finishSmokeTest(bool success, const QString &message) {
    if (!smokeTest_ || smokeFinished_) return;
    smokeFinished_ = true;
    appStatus_->setText(success ? QStringLiteral("Graphical smoke test passed")
                                : QStringLiteral("Graphical smoke test failed"));
    progressLabel_->setText(message);
    const QString screenshot = qgetenv("QCOMPILE_LAB_SCREENSHOT").isEmpty()
        ? QDir::tempPath() + QStringLiteral("/qcompile-lab-smoke.png")
        : QString::fromLocal8Bit(qgetenv("QCOMPILE_LAB_SCREENSHOT"));
    if (!QFileInfo::exists(screenshot)) {
      // Ensure a full desktop-sized render if running on headless offscreen platform
      if (width() < 1200 || height() < 800) {
        resize(1440, 920);
      }
      if (!grab().save(screenshot)) success = false;
    } else {
      progressLabel_->setText(QStringLiteral("Screenshot path already exists; not overwritten: ") + screenshot);
      success = false;
    }
    if (!success) appStatus_->setText(QStringLiteral("Graphical smoke test failed"));
    QTimer::singleShot(250, [success, message] {
      if (!success) std::cerr << "QCompile lab smoke test failed: " << message.toStdString() << std::endl;
      else std::cout << "QCompile lab smoke test passed: " << message.toStdString() << std::endl;
      QApplication::exit(success ? 0 : 1);
    });
  }

  QString root_;
  bool smokeTest_ = false;
  bool smokeFinished_ = false;
  int smokeStage_ = 0;
  std::uint64_t requestId_ = 0;
  bool jobRunning_ = false;
  std::shared_ptr<JobToken> token_;
  std::thread worker_;
  std::mutex messagesMutex_;
  std::deque<UiMessage> messages_;
  std::optional<EnergyResult> latestResult_;
  std::optional<EnergyResult> smokeFirstResult_;
  std::vector<EnergyResult> plotResults_;
  std::vector<EnergyResult> currentSweepResults_;
  QTimer pollTimer_;
  QTimer playbackTimer_;
  MoleculeView *moleculeView_ = nullptr;
  QTableWidget *gateTable_ = nullptr;
  QTableWidget *stateTable_ = nullptr;
  QLabel *stateStageLabel_ = nullptr;
  QLabel *stateNormLabel_ = nullptr;
  QSlider *stepSlider_ = nullptr;
  QPushButton *resetButton_ = nullptr;
  QPushButton *stepButton_ = nullptr;
  QPushButton *playButton_ = nullptr;
  QPushButton *pauseButton_ = nullptr;
  QDoubleSpinBox *distanceSpin_ = nullptr;
  QLabel *previewLabel_ = nullptr;
  QLabel *appStatus_ = nullptr;
  QLabel *progressLabel_ = nullptr;
  QProgressBar *progress_ = nullptr;
  QPushButton *runButton_ = nullptr;
  QPushButton *sweepButton_ = nullptr;
  QPushButton *cancelButton_ = nullptr;
  QLabel *energyGeometry_ = nullptr;
  QLabel *electronicValue_ = nullptr;
  QLabel *nuclearValue_ = nullptr;
  QLabel *totalValue_ = nullptr;
  QLabel *exactValue_ = nullptr;
  QLabel *referenceValue_ = nullptr;
  QLabel *errorValue_ = nullptr;
  QLabel *optimizerValue_ = nullptr;
  QButtonGroup *predictionGroup_ = nullptr;
  QRadioButton *predictionLower_ = nullptr;
  QRadioButton *predictionLonger_ = nullptr;
  QRadioButton *predictionEqual_ = nullptr;
  QLabel *predictionFeedback_ = nullptr;
  QLabel *explanation_ = nullptr;
  EnergyPlot *energyPlot_ = nullptr;
  QTabWidget *rightTabs_ = nullptr;
};

} // namespace

int main(int argc, char *argv[]) {
  QApplication application(argc, argv);
  const QStringList arguments = application.arguments();
  const bool smokeTest = arguments.contains(QStringLiteral("--smoke-test"));
  try {
    MainWindow window(smokeTest);
    window.showMaximized();
    return application.exec();
  } catch (const std::exception &error) {
    qCritical().noquote() << error.what();
    return 1;
  }
}
