#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"
#include "chemistry/vqe_convergence.h"
#include "simulator/simulator.h"
#include "simulator/simulation_session.h"
#include "compiler/compiler_facade.h"
#include "gui/circuit_view.h"
#include "experiment/experiment_bundle.h"
#include "compiler/circuit_dag.h"
#include <QAbstractItemView>
#include <QFileDialog>
#include <QMessageBox>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QFile>
#include <QPlainTextEdit>
#include <QTextBlock>
#include <QTextCursor>
#include <QButtonGroup>
#include <QDateTime>
#include <QDebug>
#include <QDir>
#include <QDoubleSpinBox>
#include <QEventLoop>
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
#include <QScrollBar>
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
#include <new>
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
  enum class Type { Progress, Point, Complete, Failed, Cancelled, Convergence } type;
  std::uint64_t id = 0;
  int completed = 0;
  int total = 0;
  bool sweep = false;
  QString text;
  EnergyResult result;
  VqeConvergencePoint convergence;
  QString molecule;
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
  const bool useEngine = (qgetenv("QCOMPILE_LAB_USE_ENGINE") == "1") ||
      (!QFileInfo(python).isFile() || !QFileInfo(python).isExecutable());

  if (useEngine) {
    if (!QFileInfo(benchmark).isFile() || !QFileInfo(benchmark).isExecutable()) {
      throw std::runtime_error(
          QString("qcompile_benchmarks executable is missing or not executable: %1 (build it with make qcompile_benchmarks or set QCOMPILE_LAB_BENCHMARKS to override)")
              .arg(benchmark).toStdString());
    }
    progress.publish(QString("Computing STO-3G integrals & VQE via C++ engine at %1 Angstrom…")
                         .arg(distance, 0, 'f', 3));
    const ProcessResult calculation = runProcess(
        QStringLiteral("QCompile benchmark"), benchmark,
        {QStringLiteral("--engine-h2-json"), QString::number(distance, 'g', 17)}, root, token);
    return parseBenchmarkJson(calculation.standardOutput, distance);
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

struct MoleculeSeries {
  QString name;
  QString formula;
  QColor pointColor;
  QColor exactColor;
  std::vector<EnergyResult> points;
};

class EnergyPlot final : public QWidget {
public:
  explicit EnergyPlot(QWidget *parent = nullptr) : QWidget(parent) {
    setMinimumHeight(240);
  }

  void setResults(const std::vector<EnergyResult> &results) {
    results_ = results;
    if (seriesList_.empty()) {
      seriesList_.push_back({
          QStringLiteral("H₂ (neutral singlet, STO-3G)"),
          QStringLiteral("H2"),
          QColor(56, 189, 248),
          QColor(168, 85, 247),
          results
      });
    } else {
      seriesList_[0].points = results;
    }
    update();
  }

  void setSeriesList(const std::vector<MoleculeSeries> &seriesList) {
    seriesList_ = seriesList;
    update();
  }

  void setViewMode(int mode) {
    viewMode_ = mode; // 0: Overlay (All), 1: H2 only, 2: HeH+ only
    update();
  }

  int viewMode() const { return viewMode_; }
  const std::vector<MoleculeSeries> &seriesList() const { return seriesList_; }

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

    // Title
    QString title = QStringLiteral("BOND-LENGTH ENERGY SCAN (STO-3G) // MULTI-MOLECULE PES OVERLAY (H₂ & HeH⁺)");
    if (viewMode_ == 1) title = QStringLiteral("H₂ BOND-LENGTH ENERGY SCAN (STO-3G) // DISCRETE RUN SAMPLES");
    if (viewMode_ == 2) title = QStringLiteral("HeH⁺ BENCHMARK DATASET (STO-3G) // 1.463283 BOHR (0.774 Å)");

    painter.setPen(QColor(90, 110, 133));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));
    painter.drawText(QRect(12, 10, width() - 24, 16), Qt::AlignLeft, title);

    // Collect active points across visible series
    std::vector<const MoleculeSeries *> activeSeries;
    if (viewMode_ == 0) {
      for (const auto &s : seriesList_) if (!s.points.empty()) activeSeries.push_back(&s);
    } else if (viewMode_ == 1 && !seriesList_.empty() && !seriesList_[0].points.empty()) {
      activeSeries.push_back(&seriesList_[0]);
    } else if (viewMode_ == 2 && seriesList_.size() > 1 && !seriesList_[1].points.empty()) {
      activeSeries.push_back(&seriesList_[1]);
    }

    if (activeSeries.empty()) {
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 9));
      painter.drawText(rect(), Qt::AlignCenter,
                       QStringLiteral("No calculation samples. Run a geometry calculation, 5-point sweep, or overlay."));
      return;
    }

    double minX = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

    for (const auto *s : activeSeries) {
      for (const auto &pt : s->points) {
        minX = std::min(minX, pt.distanceAngstrom);
        maxX = std::max(maxX, pt.distanceAngstrom);
        minY = std::min({minY, pt.totalEnergy, pt.exactTotalEnergy});
        maxY = std::max({maxY, pt.totalEnergy, pt.exactTotalEnergy});
      }
    }
    if (maxX - minX < 1e-4) { minX -= 0.15; maxX += 0.15; }
    if (maxY - minY < 1e-6) { minY -= 0.02; maxY += 0.02; }
    const double yPad = (maxY - minY) * 0.16;
    minY -= yPad;
    maxY += yPad;
    const double xPad = (maxX - minX) * 0.08;
    minX -= xPad;
    maxX += xPad;

    const QRectF chart(72, 34, std::max(10, width() - 92), std::max(20, height() - 95));
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

    // Draw series
    for (const auto *s : activeSeries) {
      const bool isHeH = (s->formula == QStringLiteral("HeH+"));

      // Connect points if more than 1 point
      std::vector<EnergyResult> sorted = s->points;
      std::sort(sorted.begin(), sorted.end(),
                [](const EnergyResult &a, const EnergyResult &b) {
                  return a.distanceAngstrom < b.distanceAngstrom;
                });
      if (sorted.size() > 1) {
        QColor lineCol = s->pointColor;
        lineCol.setAlpha(120);
        painter.setPen(QPen(lineCol, 1.2, Qt::SolidLine));
        for (std::size_t i = 1; i < sorted.size(); ++i) {
          painter.drawLine(point(sorted[i - 1].distanceAngstrom, sorted[i - 1].totalEnergy),
                           point(sorted[i].distanceAngstrom, sorted[i].totalEnergy));
        }
      }

      for (const auto &res : s->points) {
        const QPointF calcPt = point(res.distanceAngstrom, res.totalEnergy);
        const QPointF exactPt = point(res.distanceAngstrom, res.exactTotalEnergy);

        if (!isHeH) {
          // H2: Circular points and triangular exact marks
          painter.setPen(QPen(s->pointColor, 1.5));
          painter.setBrush(s->pointColor.darker(140));
          painter.drawEllipse(calcPt, 4.5, 4.5);

          painter.setPen(QPen(s->exactColor, 1.5));
          painter.setBrush(Qt::NoBrush);
          QPolygonF tri;
          tri << QPointF(exactPt.x(), exactPt.y() - 4.5)
              << QPointF(exactPt.x() - 4.5, exactPt.y() + 3.5)
              << QPointF(exactPt.x() + 4.5, exactPt.y() + 3.5);
          painter.drawPolygon(tri);
        } else {
          // HeH+: Square points and diamond exact marks
          painter.setPen(QPen(s->pointColor, 1.5));
          painter.setBrush(s->pointColor.darker(140));
          painter.drawRect(QRectF(calcPt.x() - 4.0, calcPt.y() - 4.0, 8.0, 8.0));

          painter.setPen(QPen(s->exactColor, 1.5));
          painter.setBrush(Qt::NoBrush);
          QPolygonF dia;
          dia << QPointF(exactPt.x(), exactPt.y() - 5.0)
              << QPointF(exactPt.x() + 5.0, exactPt.y())
              << QPointF(exactPt.x(), exactPt.y() + 5.0)
              << QPointF(exactPt.x() - 5.0, exactPt.y());
          painter.drawPolygon(dia);
        }
      }
    }

    // Axis title
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8));
    painter.drawText(QRectF(chart.left(), chart.bottom() + 18, chart.width(), 16),
                     Qt::AlignCenter, QStringLiteral("BOND LENGTH R (Å)"));

    // Legend
    painter.setPen(QColor(56, 189, 248));
    painter.drawText(QRect(12, height() - 20, 115, 16), Qt::AlignLeft,
                     QStringLiteral("● H₂ VQE"));
    painter.setPen(QColor(168, 85, 247));
    painter.drawText(QRect(130, height() - 20, 130, 16), Qt::AlignLeft,
                     QStringLiteral("△ H₂ Exact"));
    painter.setPen(QColor(16, 185, 129));
    painter.drawText(QRect(265, height() - 20, 150, 16), Qt::AlignLeft,
                     QStringLiteral("■ HeH⁺ VQE (0.774 Å)"));
    painter.setPen(QColor(245, 158, 11));
    painter.drawText(QRect(420, height() - 20, 155, 16), Qt::AlignLeft,
                     QStringLiteral("◆ HeH⁺ Exact"));
  }

private:
  std::vector<EnergyResult> results_;
  std::vector<MoleculeSeries> seriesList_;
  int viewMode_ = 0;
};

class ConvergencePlot final : public QWidget {
public:
  explicit ConvergencePlot(QWidget *parent = nullptr) : QWidget(parent) {
    setMinimumHeight(240);
  }

  void reset(double exactReference = 0.0, const QString &molecule = QStringLiteral("H2")) {
    points_.clear();
    exactReference_ = exactReference;
    molecule_ = molecule;
    statusText_ = QStringLiteral("Ready. Click 'Run Live VQE' to stream optimizer trajectory.");
    update();
  }

  void addPoint(const VqeConvergencePoint &pt) {
    points_.push_back(pt);
    update();
  }

  void setStatus(const QString &status) {
    statusText_ = status;
    update();
  }

  const std::vector<VqeConvergencePoint> &points() const { return points_; }
  double exactReference() const { return exactReference_; }

protected:
  void paintEvent(QPaintEvent *) override {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    painter.fillRect(rect(), QColor(11, 15, 23));
    painter.setPen(QPen(QColor(30, 41, 59), 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(rect().adjusted(0, 0, -1, -1));

    painter.setPen(QColor(90, 110, 133));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));
    painter.drawText(QRect(12, 10, width() - 24, 16), Qt::AlignLeft,
                     QStringLiteral("%1 LIVE VQE CONVERGENCE TRAJECTORY // REAL-TIME EVALUATION STREAM")
                         .arg(molecule_.toUpper()));

    if (points_.empty()) {
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 9));
      painter.drawText(rect(), Qt::AlignCenter, statusText_);
      return;
    }

    int minX = 1;
    int maxX = std::max(2, static_cast<int>(points_.size()));
    double minY = points_.front().totalEnergy;
    double maxY = minY;
    if (std::isfinite(exactReference_) && exactReference_ != 0.0) {
      minY = std::min(minY, exactReference_);
      maxY = std::max(maxY, exactReference_);
    }
    for (const auto &pt : points_) {
      minY = std::min(minY, pt.totalEnergy);
      maxY = std::max(maxY, pt.totalEnergy);
    }
    if (maxY - minY < 1e-8) { minY -= 0.01; maxY += 0.01; }
    const double yPad = (maxY - minY) * 0.15;
    minY -= yPad;
    maxY += yPad;

    const QRectF chart(72, 34, std::max(10, width() - 92), std::max(20, height() - 95));
    painter.fillRect(chart, QColor(15, 23, 42, 160));
    painter.setPen(QPen(QColor(42, 58, 82), 1));
    painter.drawRect(chart);

    // Horizontal grid & Y labels
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

    // Vertical grid & X labels
    for (int i = 0; i <= 4; ++i) {
      const double fraction = i / 4.0;
      const double xValue = minX + fraction * (maxX - minX);
      const double xPos = chart.left() + fraction * chart.width();
      painter.setPen(QPen(QColor(30, 45, 68), 1, Qt::DotLine));
      painter.drawLine(QPointF(xPos, chart.top()), QPointF(xPos, chart.bottom()));
      painter.setPen(QColor(100, 116, 139));
      painter.setFont(QFont(QStringLiteral("Monospace"), 8));
      painter.drawText(QRectF(xPos - 25, chart.bottom() + 4, 50, 14),
                       Qt::AlignCenter, QString::number(static_cast<int>(xValue)));
    }

    // Exact reference dashed line
    if (std::isfinite(exactReference_) && exactReference_ != 0.0) {
      const double refY = chart.bottom() - (exactReference_ - minY) / (maxY - minY) * chart.height();
      if (refY >= chart.top() && refY <= chart.bottom()) {
        painter.setPen(QPen(QColor(192, 132, 252, 200), 1.5, Qt::DashLine));
        painter.drawLine(QPointF(chart.left(), refY), QPointF(chart.right(), refY));
      }
    }

    // Trial evaluation dots
    for (const auto &pt : points_) {
      const QPointF p = point(pt.evaluation, pt.totalEnergy);
      painter.setPen(Qt::NoPen);
      painter.setBrush(QColor(56, 189, 248, 140));
      painter.drawEllipse(p, 2.5, 2.5);
    }

    // Best-so-far trajectory line
    std::vector<QPointF> trajectory;
    double currentBest = std::numeric_limits<double>::infinity();
    for (const auto &pt : points_) {
      if (pt.totalEnergy < currentBest) {
        currentBest = pt.totalEnergy;
      }
      trajectory.push_back(point(pt.evaluation, currentBest));
    }
    if (trajectory.size() > 1) {
      painter.setPen(QPen(QColor(245, 158, 11), 2, Qt::SolidLine));
      for (std::size_t i = 1; i < trajectory.size(); ++i) {
        painter.drawLine(trajectory[i - 1], trajectory[i]);
      }
    }

    // Highlight final minimum
    if (!trajectory.empty()) {
      painter.setPen(QPen(QColor(245, 158, 11), 2));
      painter.setBrush(QColor(251, 191, 36));
      painter.drawEllipse(trajectory.back(), 5, 5);
    }

    // Chart footer
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont(QStringLiteral("Monospace"), 8));
    painter.drawText(QRectF(chart.left(), chart.bottom() + 18, chart.width(), 16),
                     Qt::AlignCenter, QStringLiteral("OPTIMIZATION EVALUATION (N)"));

    painter.setPen(QColor(56, 189, 248));
    painter.drawText(QRect(12, height() - 20, 110, 16), Qt::AlignLeft,
                     QStringLiteral("● Trial Evals"));
    painter.setPen(QColor(245, 158, 11));
    painter.drawText(QRect(125, height() - 20, 140, 16), Qt::AlignLeft,
                     QStringLiteral("― Best Trajectory"));
    painter.setPen(QColor(192, 132, 252));
    painter.drawText(QRect(270, height() - 20, 130, 16), Qt::AlignLeft,
                     QStringLiteral("-- Exact Target"));

    if (!points_.empty()) {
      painter.setPen(QColor(245, 158, 11));
      painter.drawText(QRect(410, height() - 20, width() - 420, 16), Qt::AlignLeft,
                       QStringLiteral("Best: %1 Ha (Eval %2)")
                           .arg(currentBest, 0, 'f', 6)
                           .arg(points_.back().evaluation));
    }
  }

private:
  std::vector<VqeConvergencePoint> points_;
  double exactReference_ = 0.0;
  QString molecule_ = QStringLiteral("H2");
  QString statusText_ = QStringLiteral("Ready.");
};

struct ComparisonEntry {
  ExperimentBundle bundle;
  QString fileIdentity;
};

class MainWindow final : public QMainWindow {
public:
  explicit MainWindow(bool smokeTest) : root_(repositoryRoot()), smokeTest_(smokeTest) {
    setWindowTitle(QStringLiteral("QCompile — Interactive Quantum Programming & Simulation Workbench"));
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
        "QPlainTextEdit { background:#090d16; border:1px solid #1e293b; border-radius:2px; color:#f1f5f9; font-family:monospace; font-size:11px; padding:6px; selection-background-color:#1e3a5f; }"
        "QComboBox { background:#0f172a; border:1px solid #29384e; border-radius:2px; padding:3px 8px; color:#f1f5f9; font-size:11px; font-weight:600; min-width:180px; }"
        "QComboBox::drop-down { border:0px; width:18px; }"
        "QComboBox QAbstractItemView { background:#0f172a; border:1px solid #1e293b; selection-background-color:#1e3a5f; selection-color:#f1f5f9; color:#cbd5e1; }"
        "QCheckBox { spacing:6px; color:#94a3b8; font-size:11px; font-weight:600; }"
        "QCheckBox::indicator { width:13px; height:13px; border:1px solid #334155; border-radius:2px; background:#0f172a; }"
        "QCheckBox::indicator:checked { background:#0284c7; border-color:#38bdf8; }"
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
    updateEnergyPlotSeries();
    loadPresetExample(0);
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

    auto *modLabel = new QLabel(QStringLiteral("INTERACTIVE WORKBENCH"));
    modLabel->setStyleSheet(QStringLiteral("font-family:monospace;font-size:12px;font-weight:700;letter-spacing:1px;color:#f1f5f9"));

    auto *sep2 = new QLabel(QStringLiteral("·"));
    sep2->setStyleSheet(QStringLiteral("color:#334155"));

    auto *subLabel = new QLabel(QStringLiteral("Quantum Circuit Simulation & Molecular PES"));
    subLabel->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b"));

    header->addWidget(brandLabel);
    header->addWidget(sep1);
    header->addWidget(modLabel);
    header->addWidget(sep2);
    header->addWidget(subLabel);
    header->addStretch(1);

    appStatus_ = new QLabel(QStringLiteral("READY"));
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

    mainTabs_ = new QTabWidget;
    mainTabs_->addTab(buildCircuitWorkbench(), QStringLiteral("⚛  Quantum Circuit Workbench"));
    mainTabs_->addTab(buildChemistryWorkbench(), QStringLiteral("🔬  Molecular PES Laboratory"));
    mainTabs_->addTab(buildComparisonInspector(), QStringLiteral("📊  Comparison Inspector"));
    outer->addWidget(mainTabs_, 1);

    setCentralWidget(central);
  }

  QWidget *buildCircuitWorkbench() {
    auto *circuitContainer = new QWidget;
    auto *outerLayout = new QVBoxLayout(circuitContainer);
    outerLayout->setContentsMargins(4, 4, 4, 4);
    outerLayout->setSpacing(6);

    auto *splitter = new QSplitter(Qt::Horizontal);
    splitter->setChildrenCollapsible(false);

    // LEFT PANE: DSL Source Code & Diagnostics
    auto *leftBox = new QGroupBox(QStringLiteral("QUANTUM DSL EDITOR & DIAGNOSTICS"));
    auto *leftLayout = new QVBoxLayout(leftBox);
    leftLayout->setContentsMargins(8, 12, 8, 8);
    leftLayout->setSpacing(6);

    // Keep the frequently used compile action separate from secondary actions.
    auto *toolbar = new QHBoxLayout;
    toolbar->setSpacing(8);

    auto *presetLabel = new QLabel(QStringLiteral("Preset:"));
    presetLabel->setStyleSheet(QStringLiteral("font-weight:700;color:#94a3b8"));
    exampleCombo_ = new QComboBox;
    exampleCombo_->addItem(QStringLiteral("Bell State (Entanglement)"));
    exampleCombo_->addItem(QStringLiteral("Superposition & Phase (H + RZ)"));
    exampleCombo_->addItem(QStringLiteral("3-Qubit GHZ State"));
    exampleCombo_->addItem(QStringLiteral("Single Qubit Bit-Flip (X)"));
    exampleCombo_->addItem(QStringLiteral("HeH⁺ 4-Qubit Singlet VQE Circuit"));
    exampleCombo_->addItem(QStringLiteral("Syntax Error Diagnostic Demo"));

    compileBtn_ = new QPushButton(QStringLiteral("Compile Circuit"));
    compileBtn_->setObjectName(QStringLiteral("primary"));
    compileBtn_->setFixedHeight(26);

    optimizeCheck_ = new QCheckBox(QStringLiteral("Optimize Pass"));
    optimizeCheck_->setChecked(true);

    saveExpBtn_ = new QPushButton(QStringLiteral("💾 Save"));
    saveExpBtn_->setFixedHeight(26);
    openExpBtn_ = new QPushButton(QStringLiteral("📂 Open"));
    openExpBtn_->setFixedHeight(26);
    exportSvgBtn_ = new QPushButton(QStringLiteral("⤓ SVG"));
    exportSvgBtn_->setFixedHeight(26);
    exportPngBtn_ = new QPushButton(QStringLiteral("⤓ PNG"));
    exportPngBtn_->setFixedHeight(26);

    toolbar->addWidget(presetLabel);
    toolbar->addWidget(exampleCombo_, 1);
    toolbar->addWidget(compileBtn_);
    leftLayout->addLayout(toolbar);

    auto *secondaryToolbar = new QHBoxLayout;
    secondaryToolbar->setSpacing(8);
    secondaryToolbar->addWidget(optimizeCheck_);
    secondaryToolbar->addWidget(saveExpBtn_);
    secondaryToolbar->addWidget(openExpBtn_);
    secondaryToolbar->addWidget(exportSvgBtn_);
    secondaryToolbar->addWidget(exportPngBtn_);
    secondaryToolbar->addStretch(1);
    leftLayout->addLayout(secondaryToolbar);

    // Status label
    compileStatusLabel_ = new QLabel(QStringLiteral("READY"));
    compileStatusLabel_->setStyleSheet(QStringLiteral(
        "font-family:monospace;font-size:11px;font-weight:700;color:#38bdf8;padding:4px 8px;background:#0f172a;border:1px solid #1e293b;border-radius:2px"));
    leftLayout->addWidget(compileStatusLabel_);

    // Code Editor
    dslEditor_ = new QPlainTextEdit;
    dslEditor_->setLineWrapMode(QPlainTextEdit::NoWrap);
    dslEditor_->setPlaceholderText(QStringLiteral("// Enter Quantum DSL source code here..."));
    leftLayout->addWidget(dslEditor_, 3);

    // Diagnostics Table
    diagnosticsTable_ = new QTableWidget(0, 4);
    diagnosticsTable_->setHorizontalHeaderLabels({
        QStringLiteral("Line"), QStringLiteral("Column"),
        QStringLiteral("Severity"), QStringLiteral("Diagnostic Message")});
    diagnosticsTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    diagnosticsTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    diagnosticsTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    diagnosticsTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    diagnosticsTable_->verticalHeader()->hide();
    diagnosticsTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    diagnosticsTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    diagnosticsTable_->setSelectionMode(QAbstractItemView::SingleSelection);
    diagnosticsTable_->setFixedHeight(120);
    diagnosticsTable_->hide();
    leftLayout->addWidget(diagnosticsTable_, 1);

    splitter->addWidget(leftBox);

    // RIGHT PANE: Interactive Simulation & State Inspector
    auto *rightPane = new QWidget;
    auto *rightLayout = new QVBoxLayout(rightPane);
    rightLayout->setContentsMargins(0, 0, 0, 0);
    rightLayout->setSpacing(6);

    auto *circuitBox = new QGroupBox(QStringLiteral(
        "CIRCUIT VISUALIZATION // LOGICAL LAYERS, SEQUENTIAL EXECUTION"));
    auto *circuitBoxLayout = new QVBoxLayout(circuitBox);
    circuitBoxLayout->setContentsMargins(8, 12, 8, 8);
    circuitBoxLayout->setSpacing(6);

    // Control bar
    auto *controls = new QHBoxLayout;
    controls->setSpacing(6);

    circuitResetBtn_ = new QPushButton(QStringLiteral("Reset"));
    circuitStepBackBtn_ = new QPushButton(QStringLiteral("◀ Back"));
    circuitStepFwdBtn_ = new QPushButton(QStringLiteral("Next ▶"));
    circuitPlayBtn_ = new QPushButton(QStringLiteral("▶ Play"));
    circuitPauseBtn_ = new QPushButton(QStringLiteral("⏸ Pause"));
    circuitResetBtn_->setToolTip(QStringLiteral("Reset to the initial state"));
    circuitStepBackBtn_->setToolTip(QStringLiteral("Step backward"));
    circuitStepFwdBtn_->setToolTip(QStringLiteral("Step forward"));
    circuitPlayBtn_->setToolTip(QStringLiteral("Play sequential instructions"));
    circuitPauseBtn_->setToolTip(QStringLiteral("Pause playback"));
    circuitResetBtn_->setFixedHeight(24);
    circuitStepBackBtn_->setFixedHeight(24);
    circuitStepFwdBtn_->setFixedHeight(24);
    circuitPlayBtn_->setFixedHeight(24);
    circuitPauseBtn_->setFixedHeight(24);

    circuitStepSlider_ = new QSlider(Qt::Horizontal);
    circuitStepSlider_->setRange(0, 0);
    circuitStepSlider_->setFixedHeight(18);

    circuitStepLabel_ = new QLabel(QStringLiteral("STEP 0 / 0"));
    circuitStepLabel_->setStyleSheet(QStringLiteral(
        "font-family:monospace;font-size:10px;font-weight:700;color:#38bdf8;padding:2px 6px;background:#0f172a;border:1px solid #1e293b;border-radius:2px"));

    viewOptimizedRadio_ = new QRadioButton(QStringLiteral("Optimized"));
    viewRawRadio_ = new QRadioButton(QStringLiteral("Raw"));
    viewOptimizedRadio_->setChecked(true);
    auto *viewGroup = new QButtonGroup(this);
    viewGroup->addButton(viewOptimizedRadio_);
    viewGroup->addButton(viewRawRadio_);

    controls->addWidget(circuitResetBtn_);
    controls->addWidget(circuitStepBackBtn_);
    controls->addWidget(circuitStepFwdBtn_);
    controls->addWidget(circuitPlayBtn_);
    controls->addWidget(circuitPauseBtn_);
    controls->addWidget(circuitStepSlider_, 1);
    controls->addWidget(circuitStepLabel_);
    auto *viewControls = new QHBoxLayout;
    viewControls->setSpacing(6);
    viewControls->addWidget(viewOptimizedRadio_);
    viewControls->addWidget(viewRawRadio_);
    viewControls->addStretch(1);
    circuitBoxLayout->addLayout(controls);
    circuitBoxLayout->addLayout(viewControls);

    circuitView_ = new CircuitView;
    circuitScrollArea_ = new QScrollArea;
    circuitScrollArea_->setWidget(circuitView_);
    circuitScrollArea_->setWidgetResizable(true);
    circuitScrollArea_->setMinimumHeight(180);
    circuitScrollArea_->setStyleSheet(QStringLiteral("QScrollArea { border:1px solid #1e293b; background:#0b0f17; }"));
    circuitScrollArea_->viewport()->installEventFilter(this);
    circuitViewportHint_ = new QLabel;
    circuitViewportHint_->setWordWrap(true);
    circuitViewportHint_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:9px;color:#fbbf24"));
    circuitViewportHint_->hide();
    circuitBoxLayout->addWidget(circuitScrollArea_);
    circuitBoxLayout->addWidget(circuitViewportHint_);
    // Lower: Statevector & Measurement Inspector
    auto *stateBox = new QGroupBox(QStringLiteral("QUANTUM STATE & MEASUREMENT INSPECTOR"));
    auto *stateBoxLayout = new QVBoxLayout(stateBox);
    stateBoxLayout->setContentsMargins(8, 12, 8, 8);
    stateBoxLayout->setSpacing(6);

    auto *readoutRow = new QHBoxLayout;
    circuitStageLabel_ = new QLabel(QStringLiteral("Selected Stage: Initial |0...0⟩ register"));
    circuitStageLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:11px;color:#38bdf8;font-weight:600"));
    circuitMetricsLabel_ = new QLabel(QStringLiteral("Register: 0 Qubits · Norm: 1.000000"));
    circuitMetricsLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#94a3b8"));
    readoutRow->addWidget(circuitStageLabel_);
    readoutRow->addStretch(1);
    readoutRow->addWidget(circuitMetricsLabel_);
    stateBoxLayout->addLayout(readoutRow);

    circuitMeasurementBanner_ = new QLabel;
    circuitMeasurementBanner_->setWordWrap(true);
    circuitMeasurementBanner_->setStyleSheet(QStringLiteral(
        "font-family:monospace;font-size:11px;font-weight:700;color:#c7d2fe;background:#1e1b4b;border:1px solid #6366f1;padding:6px 10px;border-radius:3px"));
    circuitMeasurementBanner_->hide();
    stateBoxLayout->addWidget(circuitMeasurementBanner_);

    circuitStateTable_ = new QTableWidget(0, 4);
    circuitStateTable_->setHorizontalHeaderLabels({
        QStringLiteral("Basis |b⟩"),
        QStringLiteral("Amplitude Re + i Im"),
        QStringLiteral("Probability |c|²"),
        QStringLiteral("Visual Distribution")
    });
    circuitStateTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    circuitStateTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    circuitStateTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    circuitStateTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    circuitStateTable_->verticalHeader()->hide();
    circuitStateTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    circuitStateTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    circuitStateTable_->setSelectionMode(QAbstractItemView::SingleSelection);
    stateBoxLayout->addWidget(circuitStateTable_, 1);

    auto *circuitStateSplitter = new QSplitter(Qt::Vertical, rightPane);
    circuitStateSplitter->setChildrenCollapsible(false);
    circuitStateSplitter->addWidget(circuitBox);
    circuitStateSplitter->addWidget(stateBox);
    circuitStateSplitter->setStretchFactor(0, 1);
    circuitStateSplitter->setStretchFactor(1, 2);
    circuitStateSplitter->setSizes({350, 500});
    rightLayout->addWidget(circuitStateSplitter, 1);
    splitter->addWidget(rightPane);

    splitter->setStretchFactor(0, 9);
    splitter->setStretchFactor(1, 13);
    outerLayout->addWidget(splitter);

    // Connect signals
    connect(compileBtn_, &QPushButton::clicked, [this] { compileDsl(); });
    connect(optimizeCheck_, &QCheckBox::toggled, [this] {
      if (!isCompiling_) compileDsl();
    });
    connect(exampleCombo_, qOverload<int>(&QComboBox::currentIndexChanged), [this](int idx) { loadPresetExample(idx); });
    connect(circuitResetBtn_, &QPushButton::clicked, [this] { circuitReset(); });
    connect(circuitStepBackBtn_, &QPushButton::clicked, [this] { circuitStepBackward(); });
    connect(circuitStepFwdBtn_, &QPushButton::clicked, [this] { circuitStepForward(); });
    connect(circuitPlayBtn_, &QPushButton::clicked, [this] { circuitPlay(); });
    connect(circuitPauseBtn_, &QPushButton::clicked, [this] { circuitPause(); });
    connect(&circuitPlayTimer_, &QTimer::timeout, [this] { circuitStepForward(); });
    connect(circuitStepSlider_, &QSlider::valueChanged, [this](int val) { circuitSeek(val); });
    connect(viewOptimizedRadio_, &QRadioButton::toggled, [this] { refreshCircuitFromCompilation(); });
    connect(diagnosticsTable_, &QTableWidget::cellClicked,
            [this](int row, int) { navigateToDiagnosticRow(row); });
    connect(saveExpBtn_, &QPushButton::clicked, [this] { saveCurrentExperiment(); });
    connect(openExpBtn_, &QPushButton::clicked, [this] { openExperimentFile(); });
    connect(exportSvgBtn_, &QPushButton::clicked, [this] { exportCircuitToSvgFile(); });
    connect(exportPngBtn_, &QPushButton::clicked, [this] { exportCircuitToPngFile(); });

    dslDebounceTimer_.setSingleShot(true);
    dslDebounceTimer_.setInterval(350);
    connect(dslEditor_, &QPlainTextEdit::textChanged, [this] {
      if (!isCompiling_) dslDebounceTimer_.start();
    });
    connect(&dslDebounceTimer_, &QTimer::timeout, [this] { compileDsl(); });

    circuitPlayTimer_.setInterval(450);

    return circuitContainer;
  }

  QWidget *buildChemistryWorkbench() {
    auto *chemContainer = new QWidget;
    auto *chemLayout = new QVBoxLayout(chemContainer);
    chemLayout->setContentsMargins(0, 0, 0, 0);
    chemLayout->setSpacing(0);

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

    auto *pesTopRow = new QHBoxLayout;
    auto *pesLabel = new QLabel(QStringLiteral("PES View Mode:"));
    pesLabel->setStyleSheet(QStringLiteral("font-weight:700;color:#94a3b8;font-size:11px"));
    pesViewCombo_ = new QComboBox;
    pesViewCombo_->addItem(QStringLiteral("All Supported Molecules (H₂ & HeH⁺ Overlay)"));
    pesViewCombo_->addItem(QStringLiteral("H₂ Potential Energy Scan (STO-3G)"));
    pesViewCombo_->addItem(QStringLiteral("HeH⁺ Benchmark Dataset (STO-3G)"));
    pesTopRow->addWidget(pesLabel);
    pesTopRow->addWidget(pesViewCombo_, 1);
    plotContainerLayout->addLayout(pesTopRow);

    energyPlot_ = new EnergyPlot;
    energyPlot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    plotContainerLayout->addWidget(energyPlot_, 1);

    explanation_ = new QLabel(QStringLiteral("Variational energy is compared with the physical-sector exact ground reference of the generated STO-3G Hamiltonian. No quantum advantage is claimed."));
    explanation_->setWordWrap(true);
    explanation_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b;padding:6px 8px;border:1px solid #1e293b;background:#090d16"));
    plotContainerLayout->addWidget(explanation_);

    rightTabs->addTab(plotContainer, QStringLiteral("Energy Scan & PES Overlay"));

    // Tab 3: Real-Time VQE Convergence Trajectory
    auto *convergenceContainer = new QWidget;
    auto *convergenceLayout = new QVBoxLayout(convergenceContainer);
    convergenceLayout->setContentsMargins(8, 8, 8, 8);
    convergenceLayout->setSpacing(6);

    convergencePlot_ = new ConvergencePlot;
    convergencePlot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    convergenceLayout->addWidget(convergencePlot_, 1);

    auto *vqeCtrlBox = new QGroupBox(QStringLiteral("LIVE VQE OPTIMIZATION CONTROLS"));
    auto *vqeCtrlLayout = new QVBoxLayout(vqeCtrlBox);
    vqeCtrlLayout->setSpacing(6);

    auto *vqeBtnRow = new QHBoxLayout;
    vqeBtnRow->setSpacing(8);
    liveH2VqeBtn_ = new QPushButton(QStringLiteral("▶ Run Live H₂ VQE"));
    liveH2VqeBtn_->setObjectName(QStringLiteral("primary"));
    liveH2VqeBtn_->setFixedHeight(28);

    liveHeHVqeBtn_ = new QPushButton(QStringLiteral("▶ Run Live HeH⁺ VQE"));
    liveHeHVqeBtn_->setFixedHeight(28);

    cancelVqeBtn_ = new QPushButton(QStringLiteral("Cancel VQE"));
    cancelVqeBtn_->setFixedHeight(28);
    cancelVqeBtn_->setEnabled(false);

    exportCsvBtn_ = new QPushButton(QStringLiteral("⤓ Export CSV"));
    exportCsvBtn_->setFixedHeight(28);

    vqeBtnRow->addWidget(liveH2VqeBtn_, 2);
    vqeBtnRow->addWidget(liveHeHVqeBtn_, 2);
    vqeBtnRow->addWidget(cancelVqeBtn_, 1);
    vqeBtnRow->addWidget(exportCsvBtn_, 1);
    vqeCtrlLayout->addLayout(vqeBtnRow);

    vqeStatusLabel_ = new QLabel(QStringLiteral("Ready. Choose H₂ or HeH⁺ to stream live optimizer evaluations."));
    vqeStatusLabel_->setStyleSheet(QStringLiteral("font-family:monospace;font-size:10px;color:#64748b;padding:4px 6px;background:#090d16;border:1px solid #1e293b"));
    vqeCtrlLayout->addWidget(vqeStatusLabel_);

    convergenceLayout->addWidget(vqeCtrlBox);

    rightTabs->addTab(convergenceContainer, QStringLiteral("VQE Convergence Trajectory"));

    connect(pesViewCombo_, qOverload<int>(&QComboBox::currentIndexChanged), [this](int idx) {
      if (energyPlot_) {
        energyPlot_->setViewMode(idx);
      }
    });
    connect(liveH2VqeBtn_, &QPushButton::clicked, [this] { startLiveVqe(false); });
    connect(liveHeHVqeBtn_, &QPushButton::clicked, [this] { startLiveVqe(true); });
    connect(cancelVqeBtn_, &QPushButton::clicked, [this] {
      if (token_) token_->cancelled.store(true);
    });
    connect(exportCsvBtn_, &QPushButton::clicked, [this] { exportConvergenceToCsvFile(); });
    columns->addWidget(rightTabs);

    columns->setStretchFactor(0, 11);
    columns->setStretchFactor(1, 10);

    chemLayout->addWidget(columns, 1);
    return chemContainer;
  }
  QWidget *buildComparisonInspector() {
    auto *compContainer = new QWidget;
    auto *layout = new QVBoxLayout(compContainer);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(8);

    // Top control box
    auto *ctrlBox = new QGroupBox(QStringLiteral("LOADED EXPERIMENT BUNDLES (.QCOMP)"));
    auto *ctrlLayout = new QVBoxLayout(ctrlBox);
    ctrlLayout->setSpacing(6);

    auto *bar = new QHBoxLayout;
    auto *addBtn = new QPushButton(QStringLiteral("＋ Add Experiment Bundle (.qcomp)"));
    addBtn->setObjectName(QStringLiteral("primary"));
    auto *clearBtn = new QPushButton(QStringLiteral("Clear List"));
    auto *infoLabel = new QLabel(QStringLiteral("Load 2 or more .qcomp bundles to compare physical properties, energy deltas, and circuit metrics."));
    infoLabel->setStyleSheet(QStringLiteral("color:#94a3b8;font-size:11px"));

    bar->addWidget(addBtn);
    bar->addWidget(clearBtn);
    bar->addSpacing(12);
    bar->addWidget(infoLabel, 1);
    ctrlLayout->addLayout(bar);

    comparisonListTable_ = new QTableWidget(0, 7);
    comparisonListTable_->setHorizontalHeaderLabels({
        QStringLiteral("Title / File Identity"),
        QStringLiteral("Species"),
        QStringLiteral("Qubits"),
        QStringLiteral("Instructions"),
        QStringLiteral("Geometry"),
        QStringLiteral("Total Energy (Ha)"),
        QStringLiteral("Evals")
    });
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(4, QHeaderView::ResizeToContents);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(5, QHeaderView::ResizeToContents);
    comparisonListTable_->horizontalHeader()->setSectionResizeMode(6, QHeaderView::ResizeToContents);
    comparisonListTable_->verticalHeader()->hide();
    comparisonListTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    comparisonListTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    comparisonListTable_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    comparisonListTable_->setFixedHeight(140);
    ctrlLayout->addWidget(comparisonListTable_);

    layout->addWidget(ctrlBox, 0);

    // Bottom comparison result box
    auto *diffBox = new QGroupBox(QStringLiteral("PAIRWISE DELTA & COMPATIBILITY ANALYSIS"));
    auto *diffLayout = new QVBoxLayout(diffBox);
    diffLayout->setSpacing(6);

    comparisonSummaryLabel_ = new QLabel(QStringLiteral("Select two experiments in the list above to view scientific comparison."));
    comparisonSummaryLabel_->setStyleSheet(QStringLiteral(
        "font-family:monospace;font-size:11px;color:#38bdf8;padding:4px 8px;background:#0f172a;border:1px solid #1e293b;border-radius:2px"));
    diffLayout->addWidget(comparisonSummaryLabel_);
    comparisonSummaryLabel_->setWordWrap(true);

    comparisonDeltaTable_ = new QTableWidget(0, 6);
    comparisonDeltaTable_->setHorizontalHeaderLabels({
        QStringLiteral("Metric / Property"),
        QStringLiteral("Experiment A"),
        QStringLiteral("Experiment B"),
        QStringLiteral("Difference (B − A)"),
        QStringLiteral("Compatibility"),
        QStringLiteral("Scientific Note")
    });
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(4, QHeaderView::ResizeToContents);
    comparisonDeltaTable_->horizontalHeader()->setSectionResizeMode(5, QHeaderView::Stretch);
    comparisonDeltaTable_->verticalHeader()->hide();
    comparisonDeltaTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    comparisonDeltaTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    diffLayout->addWidget(comparisonDeltaTable_, 1);

    layout->addWidget(diffBox, 1);

    connect(addBtn, &QPushButton::clicked, [this] {
      const QString file = QFileDialog::getOpenFileName(this, QStringLiteral("Select Experiment Bundle"), QString(), QStringLiteral("QCompile Experiment (*.qcomp);;All Files (*)"));
      if (!file.isEmpty()) addExperimentToComparison(file);
    });

    connect(clearBtn, &QPushButton::clicked, [this] {
      comparisonBundles_.clear();
      comparisonListTable_->setRowCount(0);
      comparisonDeltaTable_->setRowCount(0);
      comparisonDeltaTable_->horizontalHeaderItem(1)->setText(QStringLiteral("Experiment A"));
      comparisonDeltaTable_->horizontalHeaderItem(2)->setText(QStringLiteral("Experiment B"));
      comparisonSummaryLabel_->setText(QStringLiteral("Comparison list cleared. Add .qcomp bundles to compare."));
    });

    connect(comparisonListTable_, &QTableWidget::itemSelectionChanged, [this] {
      updateComparisonView();
    });

    return compContainer;
  }

  void clearLoadedChemistryPresentation() {
    playbackTimer_.stop();
    latestResult_.reset();
    plotResults_.clear();
    currentSweepResults_.clear();

    if (stepSlider_) stepSlider_->setValue(0);
    if (energyGeometry_) energyGeometry_->setText(QStringLiteral("—"));
    if (electronicValue_) electronicValue_->setText(QStringLiteral("—"));
    if (nuclearValue_) nuclearValue_->setText(QStringLiteral("—"));
    if (totalValue_) totalValue_->setText(QStringLiteral("—"));
    if (exactValue_) exactValue_->setText(QStringLiteral("—"));
    if (referenceValue_) referenceValue_->setText(QStringLiteral("—"));
    if (errorValue_) errorValue_->setText(QStringLiteral("—"));
    if (optimizerValue_) optimizerValue_->setText(QStringLiteral("No saved chemistry result"));
    if (progressLabel_) progressLabel_->setText(QStringLiteral("No chemistry result saved in this bundle."));
    if (predictionFeedback_) predictionFeedback_->setText(QStringLiteral("No sweep results are stored in this bundle."));
    if (explanation_) explanation_->setText(QStringLiteral("No saved sweep comparison is available."));
    if (convergencePlot_) convergencePlot_->reset(0.0, QStringLiteral("H2"));
    if (vqeStatusLabel_) vqeStatusLabel_->setText(QStringLiteral("No saved convergence trajectory."));
    if (energyPlot_) updateEnergyPlotSeries();

    refreshCircuit();
    refreshState();
    updatePreviewLabels();
  }

  void restoreSavedChemistryPresentation(const ExperimentBundle &bundle) {
    if (bundle.hasChemistry) {
      const QString formula = QString::fromStdString(bundle.formula);
      const QString units = QString::fromStdString(bundle.distanceUnits);
      constexpr double kBohrToAngstrom = 0.52917721092;
      double distanceAngstrom = 0.0;
      double distanceBohr = 0.0;
      const bool hasAngstromDistance = bundle.distanceUnits == "Angstrom";
      const bool hasBohrDistance = bundle.distanceUnits == "bohr";
      if (hasAngstromDistance) {
        distanceAngstrom = bundle.distance;
        distanceBohr = bundle.distance / kBohrToAngstrom;
      } else if (hasBohrDistance) {
        distanceBohr = bundle.distance;
        distanceAngstrom = bundle.distance * kBohrToAngstrom;
      }

      if (energyGeometry_) {
        if (hasAngstromDistance) {
          energyGeometry_->setText(QStringLiteral("%1 Å (%2 bohr)")
                                       .arg(distanceAngstrom, 0, 'f', 6)
                                       .arg(distanceBohr, 0, 'f', 9));
        } else if (hasBohrDistance) {
          energyGeometry_->setText(QStringLiteral("%1 bohr (%2 Å)")
                                       .arg(distanceBohr, 0, 'f', 9)
                                       .arg(distanceAngstrom, 0, 'f', 6));
        } else {
          energyGeometry_->setText(QStringLiteral("%1 %2")
                                       .arg(bundle.distance, 0, 'g', 12)
                                       .arg(units.isEmpty() ? QStringLiteral("(unit not saved)") : units));
        }
      }
      if (electronicValue_) electronicValue_->setText(formattedEnergy(bundle.electronicEnergy));
      if (nuclearValue_) nuclearValue_->setText(formattedEnergy(bundle.nuclearRepulsion));
      if (totalValue_) totalValue_->setText(formattedEnergy(bundle.totalEnergy));
      if (exactValue_) exactValue_->setText(formattedEnergy(bundle.exactElectronicEnergy));
      if (referenceValue_) {
        referenceValue_->setText(bundle.referenceElectronicEnergy
                                     ? formattedEnergy(*bundle.referenceElectronicEnergy)
                                     : QStringLiteral("Not saved in this bundle"));
      }
      if (errorValue_) errorValue_->setText(formattedEnergy(bundle.variationalError));
      if (optimizerValue_) {
        const QString parameter = bundle.optimizerParameterRadians
            ? QString::number(*bundle.optimizerParameterRadians, 'f', 9) + QStringLiteral(" rad")
            : QStringLiteral("not saved");
        optimizerValue_->setText(QStringLiteral("%1 evaluations · %2 · θ=%3")
                                     .arg(bundle.energyEvaluations)
                                     .arg(bundle.converged ? QStringLiteral("converged") : QStringLiteral("not converged"))
                                     .arg(parameter));
      }

      EnergyResult result;
      result.distanceAngstrom = distanceAngstrom;
      result.distanceBohr = distanceBohr;
      result.parameter = bundle.optimizerParameterRadians.value_or(0.0);
      result.electronicEnergy = bundle.electronicEnergy;
      result.exactElectronicEnergy = bundle.exactElectronicEnergy;
      result.referenceElectronicEnergy = bundle.referenceElectronicEnergy.value_or(0.0);
      result.variationalError = bundle.variationalError;
      result.nuclearRepulsionEnergy = bundle.nuclearRepulsion;
      result.totalEnergy = bundle.totalEnergy;
      result.exactTotalEnergy = bundle.exactTotalEnergy;
      result.evaluations = bundle.energyEvaluations;
      result.converged = bundle.converged;
      if (bundle.formula == "H2" && (hasAngstromDistance || hasBohrDistance)) {
        plotResults_.push_back(result);
        if (bundle.optimizerParameterRadians) latestResult_ = result;
        if (distanceSpin_) distanceSpin_->setValue(distanceAngstrom);
      }
      if (progressLabel_) {
        progressLabel_->setText(QStringLiteral("Loaded saved %1 chemistry summary; values were not recalculated.")
                                    .arg(formula.isEmpty() ? QStringLiteral("molecule") : formula));
      }
    } else if (progressLabel_) {
      progressLabel_->setText(QStringLiteral("No chemistry result was saved; no prior result is being shown."));
    }

    if (convergencePlot_) {
      const QString molecule = QString::fromStdString(
          !bundle.formula.empty() ? bundle.formula :
          (!bundle.moleculeName.empty() ? bundle.moleculeName : "H2"));
      convergencePlot_->reset(bundle.hasChemistry ? bundle.exactTotalEnergy : 0.0, molecule);
      if (bundle.hasConvergence) {
        for (const auto &point : bundle.convergenceTrajectory) convergencePlot_->addPoint(point);
      }
    }
    if (vqeStatusLabel_) {
      vqeStatusLabel_->setText(bundle.hasConvergence
                                   ? QStringLiteral("Loaded saved convergence data; optimizer was not rerun.")
                                   : QStringLiteral("No convergence trajectory was saved."));
    }
    if (energyPlot_) updateEnergyPlotSeries();
    refreshCircuit();
    refreshState();
    updatePreviewLabels();
  }

  bool saveCurrentExperiment(const QString &suggestedPath = QString()) {
    QString path = suggestedPath;
    if (path.isEmpty()) {
      path = QFileDialog::getSaveFileName(
          this, QStringLiteral("Save Experiment Bundle"),
          QDir::currentPath() + QStringLiteral("/experiment.qcomp"),
          QStringLiteral("QCompile Experiment (*.qcomp);;All Files (*)"));
    }
    if (path.isEmpty()) return false;

    try {
      ExperimentBundle b;
      b.id = QStringLiteral("exp-%1").arg(QDateTime::currentMSecsSinceEpoch()).toStdString();
      b.title = dslEditor_ ? dslEditor_->document()->firstBlock().text().trimmed().toStdString() : "Quantum Experiment";
      b.timestamp = QDateTime::currentDateTimeUtc().toString(Qt::ISODate).toStdString();
      b.description = "QCompile Qt6 workbench saved session";
      b.dslSource = dslEditor_ ? dslEditor_->toPlainText().toStdString() : "";
      b.optimize = optimizeCheck_ ? optimizeCheck_->isChecked() : true;
      for (const auto &diag : lastCompilation_.diagnostics) {
        b.diagnostics.push_back(diag.message);
      }
      b.qubitCount = circuitSession_ ? circuitSession_->qubitCount() : static_cast<int>(lastCompilation_.qubitCount);
      b.rawInstructions = lastCompilation_.rawInstructions;
      b.instructions = circuitSession_ ? circuitSession_->instructions() : lastCompilation_.instructions;

      if (circuitSession_) {
        b.currentStep = circuitSession_->currentStep();
        b.totalSteps = circuitSession_->totalSteps();
        b.isComplete = circuitSession_->isComplete();
        b.stateNorm = circuitSession_->stateNorm();
        b.statevector = circuitSession_->currentState();
        b.basisProbabilities = circuitSession_->basisProbabilities();
        for (std::size_t s = 1; s <= b.currentStep; ++s) {
          auto res = circuitSession_->measurementResultAt(s);
          if (res.has_value()) {
            b.recordedMeasurements.emplace_back(s, *res);
          }
        }
      }

      if (latestResult_) {
        b.hasChemistry = true;
        b.moleculeName = "H2";
        b.formula = "H2";
        b.basis = "STO-3G";
        b.mapping = "parity";
        b.distance = latestResult_->distanceAngstrom;
        b.distanceUnits = "Angstrom";
        b.electronicEnergy = latestResult_->electronicEnergy;
        b.nuclearRepulsion = latestResult_->nuclearRepulsionEnergy;
        b.totalEnergy = latestResult_->totalEnergy;
        b.exactElectronicEnergy = latestResult_->exactElectronicEnergy;
        b.exactTotalEnergy = latestResult_->exactTotalEnergy;
        b.variationalError = latestResult_->variationalError;
        b.energyEvaluations = latestResult_->evaluations;
        b.converged = latestResult_->converged;
        b.optimizerParameterRadians = latestResult_->parameter;
        b.referenceElectronicEnergy = latestResult_->referenceElectronicEnergy;
      }

      if (convergencePlot_ && !convergencePlot_->points().empty()) {
        b.hasConvergence = true;
        b.convergenceTrajectory = convergencePlot_->points();
      }

      saveExperimentBundle(path.toStdString(), b);
      compileStatusLabel_->setText(QStringLiteral("SAVED · %1").arg(QFileInfo(path).fileName()));
      return true;
    } catch (const std::exception &err) {
      const QString message = QString::fromLocal8Bit(err.what());
      compileStatusLabel_->setText(QStringLiteral("SAVE FAILED · ") + message);
      compileStatusLabel_->setStyleSheet(QStringLiteral("color:#f43f5e;font-weight:bold;"));
      if (appStatus_) appStatus_->setText(QStringLiteral("Experiment was not saved: ") + message);
      if (!smokeTest_) {
        QMessageBox::critical(this, QStringLiteral("Save Error"), message);
      }
      return false;
    }
  }

  bool openExperimentFile(const QString &specifiedPath = QString()) {
    QString path = specifiedPath;
    if (path.isEmpty()) {
      path = QFileDialog::getOpenFileName(
          this, QStringLiteral("Open Experiment Bundle"),
          QDir::currentPath(),
          QStringLiteral("QCompile Experiment (*.qcomp);;All Files (*)"));
    }
    if (path.isEmpty()) return false;

    bool loadedSourceApplied = false;
    try {
      const ExperimentBundle bundle = loadExperimentBundle(path.toStdString());
      loadedSourceApplied = true;
      circuitPlayTimer_.stop();
      clearLoadedChemistryPresentation();

      isCompiling_ = true;
      if (dslEditor_) dslEditor_->setPlainText(QString::fromStdString(bundle.dslSource));
      if (optimizeCheck_) optimizeCheck_->setChecked(bundle.optimize);
      isCompiling_ = false;

      compileDsl();
      if (!lastCompilation_.success) {
        if (appStatus_) {
          appStatus_->setText(QStringLiteral("Bundle source was loaded, but compilation failed; previous run output was cleared."));
        }
        return false;
      }
      if (!circuitSession_) {
        if (appStatus_) {
          appStatus_->setText(QStringLiteral("Bundle source compiled, but simulation could not initialize; previous run output was cleared."));
        }
        return false;
      }
      if (bundle.currentStep > circuitSession_->totalSteps()) {
        const QString message = QStringLiteral("saved execution step %1 exceeds the compiled circuit length %2")
                                    .arg(bundle.currentStep)
                                    .arg(circuitSession_->totalSteps());
        compileStatusLabel_->setText(QStringLiteral("OPEN FAILED · ") + message);
        if (appStatus_) appStatus_->setText(QStringLiteral("Bundle source is loaded; execution remains at step 0."));
        return false;
      }

      try {
        circuitSession_->seek(bundle.currentStep);
      } catch (const std::exception &error) {
        handleCircuitExecutionFailure(error);
        if (appStatus_) appStatus_->setText(QStringLiteral("Bundle source loaded, but replay failed; the last valid simulation step is shown."));
        return false;
      }
      updateCircuitExecutionDisplay();

      bool replayMatches = true;
      if (!bundle.statevector.empty()) {
        const auto &state = circuitSession_->currentState();
        replayMatches = bundle.totalSteps == circuitSession_->totalSteps() &&
                        bundle.isComplete == circuitSession_->isComplete() &&
                        circuitSession_->currentStep() == bundle.currentStep &&
                        state.size() == bundle.statevector.size() &&
                        std::abs(circuitSession_->stateNorm() - bundle.stateNorm) <= 1e-12;
        for (std::size_t i = 0; replayMatches && i < state.size(); ++i) {
          replayMatches = std::abs(state[i] - bundle.statevector[i]) <= 1e-12;
        }
      }
      if (replayMatches && !bundle.basisProbabilities.empty()) {
        const auto probabilities = circuitSession_->basisProbabilities();
        replayMatches = probabilities.size() == bundle.basisProbabilities.size();
        for (std::size_t i = 0; replayMatches && i < probabilities.size(); ++i) {
          replayMatches = std::abs(probabilities[i] - bundle.basisProbabilities[i]) <= 1e-12;
        }
      }
      for (const auto &measurement : bundle.recordedMeasurements) {
        if (!replayMatches) break;
        const auto replayed = circuitSession_->measurementResultAt(measurement.first);
        replayMatches = replayed.has_value() && *replayed == measurement.second;
      }
      if (!replayMatches) {
        compileStatusLabel_->setText(QStringLiteral("OPEN FAILED · saved simulation data do not match deterministic replay"));
        if (appStatus_) appStatus_->setText(QStringLiteral("Bundle source is loaded; the newly replayed state is shown, not the inconsistent saved snapshot."));
        return false;
      }

      restoreSavedChemistryPresentation(bundle);
      compileStatusLabel_->setText(QStringLiteral("OPENED · %1").arg(QFileInfo(path).fileName()));
      if (appStatus_) {
        appStatus_->setText(QStringLiteral("Saved experiment opened; source recompiled and execution replayed. Chemistry values were not recalculated."));
      }
      return true;
    } catch (const std::exception &err) {
      isCompiling_ = false;
      const QString message = QString::fromLocal8Bit(err.what());
      compileStatusLabel_->setText(QStringLiteral("OPEN FAILED · ") + message);
      if (appStatus_) {
        appStatus_->setText(loadedSourceApplied
                                ? QStringLiteral("Loaded source remains visible; unrelated previous run output was cleared.")
                                : QStringLiteral("Bundle was not opened; the existing experiment remains unchanged."));
      }
      if (!smokeTest_) {
        QMessageBox::critical(this, QStringLiteral("Open Error"), message);
      }
      return false;
    }
  }


  bool exportCircuitToSvgFile(const QString &suggestedPath = QString()) {
    QString path = suggestedPath;
    if (path.isEmpty()) {
      path = QFileDialog::getSaveFileName(
          this, QStringLiteral("Export Circuit to SVG"),
          QDir::currentPath() + QStringLiteral("/circuit.svg"),
          QStringLiteral("Scalable Vector Graphics (*.svg);;All Files (*)"));
    }
    if (path.isEmpty()) return false;

    try {
      const auto &instrs = circuitSession_ ? circuitSession_->instructions() : lastCompilation_.instructions;
      const int qCount = circuitSession_ ? circuitSession_->qubitCount() : static_cast<int>(lastCompilation_.qubitCount);
      saveCircuitSvg(path.toStdString(), instrs, qCount);
      return true;
    } catch (const std::exception &err) {
      if (!smokeTest_) QMessageBox::critical(this, QStringLiteral("Export Error"), QString::fromStdString(err.what()));
      return false;
    }
  }

  bool exportCircuitToPngFile(const QString &suggestedPath = QString()) {
    QString path = suggestedPath;
    if (path.isEmpty()) {
      path = QFileDialog::getSaveFileName(
          this, QStringLiteral("Export Circuit to PNG"),
          QDir::currentPath() + QStringLiteral("/circuit.png"),
          QStringLiteral("Portable Network Graphics (*.png);;All Files (*)"));
    }
    if (path.isEmpty()) return false;

    if (!circuitView_) return false;
    return circuitView_->exportPng(path, 2);
  }

  bool exportConvergenceToCsvFile(const QString &suggestedPath = QString()) {
    QString path = suggestedPath;
    if (path.isEmpty()) {
      path = QFileDialog::getSaveFileName(
          this, QStringLiteral("Export VQE Convergence to CSV"),
          QDir::currentPath() + QStringLiteral("/convergence.csv"),
          QStringLiteral("Comma-Separated Values (*.csv);;All Files (*)"));
    }
    if (path.isEmpty()) return false;

    try {
      ExperimentBundle b;
      b.hasChemistry = (latestResult_.has_value());
      b.nuclearRepulsion = latestResult_ ? latestResult_->nuclearRepulsionEnergy : 0.0;
      b.hasConvergence = (convergencePlot_ && !convergencePlot_->points().empty());
      if (convergencePlot_) {
        b.convergenceTrajectory = convergencePlot_->points();
      }
      saveConvergenceCsv(path.toStdString(), b);
      return true;
    } catch (const std::exception &err) {
      if (!smokeTest_) QMessageBox::critical(this, QStringLiteral("CSV Export Error"), QString::fromStdString(err.what()));
      return false;
    }
  }

  QString comparisonEntryIdentity(const ComparisonEntry &entry) const {
    if (!entry.bundle.id.empty()) {
      const auto duplicateIds = std::count_if(
          comparisonBundles_.begin(), comparisonBundles_.end(),
          [&entry](const ComparisonEntry &candidate) {
            return candidate.bundle.id == entry.bundle.id;
          });
      if (duplicateIds == 1) {
        return QStringLiteral("bundle %1").arg(QString::fromStdString(entry.bundle.id));
      }
    }
    return entry.fileIdentity;
  }

  QString comparisonEntryLabel(const ComparisonEntry &entry) const {
    const QString title = entry.bundle.title.empty()
        ? QFileInfo(entry.fileIdentity).fileName()
        : QString::fromStdString(entry.bundle.title);
    return QStringLiteral("%1 — %2").arg(title, comparisonEntryIdentity(entry));
  }

  void addExperimentToComparison(const QString &path, const std::optional<ExperimentBundle> &preloaded = std::nullopt) {
    try {
      ExperimentBundle bundle = preloaded.has_value()
          ? *preloaded : loadExperimentBundle(path.toStdString());
      const QFileInfo fileInfo(path);
      QString fileIdentity = fileInfo.canonicalFilePath();
      if (fileIdentity.isEmpty()) fileIdentity = fileInfo.absoluteFilePath();
      comparisonBundles_.push_back({std::move(bundle), fileIdentity});

      const int row = comparisonListTable_->rowCount();
      comparisonListTable_->insertRow(row);

      const ComparisonEntry &entry = comparisonBundles_.back();
      auto *identityItem = new QTableWidgetItem(comparisonEntryLabel(entry));
      identityItem->setToolTip(QStringLiteral("Bundle file: %1").arg(entry.fileIdentity));
      comparisonListTable_->setItem(row, 0, identityItem);
      comparisonListTable_->setItem(row, 1, new QTableWidgetItem(
          entry.bundle.hasChemistry ? QString::fromStdString(entry.bundle.formula) : QStringLiteral("Circuit")));
      comparisonListTable_->setItem(row, 2, new QTableWidgetItem(QString::number(entry.bundle.qubitCount)));
      comparisonListTable_->setItem(row, 3, new QTableWidgetItem(QString::number(entry.bundle.instructions.size())));
      const QString geo = entry.bundle.hasChemistry
          ? QStringLiteral("%1 %2").arg(entry.bundle.distance, 0, 'f', 3)
                .arg(QString::fromStdString(entry.bundle.distanceUnits))
          : QStringLiteral("—");
      comparisonListTable_->setItem(row, 4, new QTableWidgetItem(geo));
      const QString eStr = entry.bundle.hasChemistry
          ? QString::number(entry.bundle.totalEnergy, 'f', 6) : QStringLiteral("—");
      comparisonListTable_->setItem(row, 5, new QTableWidgetItem(eStr));
      const QString evStr = entry.bundle.hasChemistry
          ? QString::number(entry.bundle.energyEvaluations) : QStringLiteral("—");
      comparisonListTable_->setItem(row, 6, new QTableWidgetItem(evStr));

      for (std::size_t i = 0; i < comparisonBundles_.size(); ++i) {
        auto *identityItem = comparisonListTable_->item(static_cast<int>(i), 0);
        if (!identityItem) continue;
        const ComparisonEntry &loaded = comparisonBundles_[i];
        identityItem->setText(comparisonEntryLabel(loaded));
        identityItem->setToolTip(QStringLiteral("Bundle file: %1").arg(loaded.fileIdentity));
      }

      if (comparisonBundles_.size() >= 2) {
        updateComparisonView();
      }
    } catch (const std::exception &err) {
      if (!smokeTest_) QMessageBox::warning(this, QStringLiteral("Load Error"), QString::fromStdString(err.what()));
    }
  }

  void updateComparisonView() {
    if (comparisonBundles_.size() < 2) {
      comparisonDeltaTable_->setRowCount(0);
      comparisonDeltaTable_->horizontalHeaderItem(1)->setText(QStringLiteral("Experiment A"));
      comparisonDeltaTable_->horizontalHeaderItem(2)->setText(QStringLiteral("Experiment B"));
      comparisonSummaryLabel_->setText(QStringLiteral("Load at least two experiments to view delta comparisons."));
      return;
    }

    int idxA = 0;
    int idxB = 1;
    auto selectedRows = comparisonListTable_->selectionModel()->selectedRows();
    if (selectedRows.size() >= 2) {
      idxA = selectedRows[0].row();
      idxB = selectedRows[1].row();
    }

    if (idxA < 0 || idxA >= static_cast<int>(comparisonBundles_.size()) ||
        idxB < 0 || idxB >= static_cast<int>(comparisonBundles_.size())) {
      return;
    }

    const ComparisonEntry &entryA = comparisonBundles_[static_cast<std::size_t>(idxA)];
    const ComparisonEntry &entryB = comparisonBundles_[static_cast<std::size_t>(idxB)];
    const auto comp = compareExperiments(entryA.bundle, entryB.bundle);
    const QString labelA = comparisonEntryLabel(entryA);
    const QString labelB = comparisonEntryLabel(entryB);

    comparisonSummaryLabel_->setText(
        QStringLiteral("A: %1\nB: %2\n%3")
            .arg(labelA, labelB, QString::fromStdString(comp.summary)));
    comparisonDeltaTable_->horizontalHeaderItem(1)->setText(
        QStringLiteral("A · %1").arg(labelA));
    comparisonDeltaTable_->horizontalHeaderItem(2)->setText(
        QStringLiteral("B · %1").arg(labelB));
    comparisonDeltaTable_->horizontalHeaderItem(1)->setToolTip(entryA.fileIdentity);
    comparisonDeltaTable_->horizontalHeaderItem(2)->setToolTip(entryB.fileIdentity);
    comparisonDeltaTable_->setRowCount(static_cast<int>(comp.items.size()));

    for (std::size_t i = 0; i < comp.items.size(); ++i) {
      const auto &item = comp.items[i];
      const int r = static_cast<int>(i);

      auto *propItem = new QTableWidgetItem(QString::fromStdString(item.propertyName));
      auto *aItem = new QTableWidgetItem(QString::fromStdString(item.valueA));
      auto *bItem = new QTableWidgetItem(QString::fromStdString(item.valueB));
      auto *deltaItem = new QTableWidgetItem(QString::fromStdString(item.delta));
      auto *compatItem = new QTableWidgetItem(item.isCompatible ? QStringLiteral("✓ Compatible") : QStringLiteral("⚠ Incompatible"));
      compatItem->setForeground(item.isCompatible ? QColor(56, 189, 248) : QColor(244, 63, 94));
      auto *noteItem = new QTableWidgetItem(QString::fromStdString(item.note));

      comparisonDeltaTable_->setItem(r, 0, propItem);
      comparisonDeltaTable_->setItem(r, 1, aItem);
      comparisonDeltaTable_->setItem(r, 2, bItem);
      comparisonDeltaTable_->setItem(r, 3, deltaItem);
      comparisonDeltaTable_->setItem(r, 4, compatItem);
      comparisonDeltaTable_->setItem(r, 5, noteItem);
    }
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
  QString formatHeHDsl(double theta, double phi) {
    const auto instructions = hehSingletTrialCircuit(theta, phi);
    QString text = QStringLiteral(
        "# HeH+ 4-Qubit Singlet VQE Trial Circuit\n"
        "# Model: STO-3G, Jordan-Wigner mapping, 4 qubits (q0..q3)\n"
        "# State preparation: |0101⟩ followed by two singlet Givens rotations\n"
        "# Optimal parameters: theta = %1 rad, phi = %2 rad\n"
        "qubit q0, q1, q2, q3\n")
        .arg(theta, 0, 'f', 6)
        .arg(phi, 0, 'f', 6);

    for (const auto &instr : instructions) {
      if (instr.type == InstructionType::Measure) {
        text += QStringLiteral("measure");
        for (std::size_t i = 0; i < instr.qubitIndices.size(); ++i) {
          text += (i == 0 ? QStringLiteral(" q") : QStringLiteral(", q")) + QString::number(instr.qubitIndices[i]);
        }
        text += QStringLiteral("\n");
      } else if (instr.gateName == "RZ") {
        text += QStringLiteral("RZ(") + QString::number(instr.parameterValue, 'g', 17) + QStringLiteral(") q") + QString::number(instr.qubitIndices[0]) + QStringLiteral("\n");
      } else if (instr.gateName == "CNOT") {
        text += QStringLiteral("CNOT q%1, q%2\n").arg(instr.qubitIndices[0]).arg(instr.qubitIndices[1]);
      } else {
        text += QString::fromStdString(instr.gateName) + QStringLiteral(" q%1\n").arg(instr.qubitIndices[0]);
      }
    }
    text += QStringLiteral("measure q0, q1, q2, q3\n");
    return text;
  }

  EnergyResult hehBenchmarkResult() {
    EnergyResult hehRes;
    hehRes.distanceBohr = 1.463283;
    hehRes.distanceAngstrom = 1.463283 * 0.52917721092; // 0.774331 Angstrom
    hehRes.electronicEnergy = -4.218264154490334;
    hehRes.exactElectronicEnergy = -4.218264154490395;
    hehRes.nuclearRepulsionEnergy = 1.366789609392032;
    hehRes.totalEnergy = -2.851474545098302;
    hehRes.exactTotalEnergy = -2.851474545098363;
    hehRes.converged = true;
    hehRes.evaluations = 6803;
    return hehRes;
  }

  void updateEnergyPlotSeries() {
    if (!energyPlot_) return;
    MoleculeSeries h2Series;
    h2Series.name = QStringLiteral("H₂ (neutral singlet, STO-3G)");
    h2Series.formula = QStringLiteral("H2");
    h2Series.pointColor = QColor(56, 189, 248);
    h2Series.exactColor = QColor(168, 85, 247);
    h2Series.points = plotResults_;

    MoleculeSeries hehSeries;
    hehSeries.name = QStringLiteral("HeH⁺ (cation singlet, STO-3G)");
    hehSeries.formula = QStringLiteral("HeH+");
    hehSeries.pointColor = QColor(16, 185, 129);
    hehSeries.exactColor = QColor(245, 158, 11);
    hehSeries.points = { hehBenchmarkResult() };

    energyPlot_->setSeriesList({ h2Series, hehSeries });
  }

  void startLiveVqe(bool isHeH) {
    if (jobRunning_) return;
    if (worker_.joinable()) worker_.join();

    jobRunning_ = true;
    token_ = std::make_shared<JobToken>();
    token_->id = ++requestId_;
    token_->cancelled.store(false);

    if (liveH2VqeBtn_) liveH2VqeBtn_->setEnabled(false);
    if (liveHeHVqeBtn_) liveHeHVqeBtn_->setEnabled(false);
    if (cancelVqeBtn_) cancelVqeBtn_->setEnabled(true);
    if (runButton_) runButton_->setEnabled(false);
    if (sweepButton_) sweepButton_->setEnabled(false);

    const double exactRef = isHeH ? HEH_TOTAL_REFERENCE_HARTREE : H2_TOTAL_REFERENCE_HARTREE;
    const QString mol = isHeH ? QStringLiteral("HeH+") : QStringLiteral("H2");
    if (convergencePlot_) {
      convergencePlot_->reset(exactRef, mol);
    }
    if (vqeStatusLabel_) {
      vqeStatusLabel_->setText(QStringLiteral("Optimizing %1 VQE... streaming real-time evaluations").arg(mol));
    }
    if (rightTabs_) {
      rightTabs_->setCurrentIndex(2);
    }

    auto token = token_;
    worker_ = std::thread([this, isHeH, token] {
      try {
        VqeCallback callback = [this, isHeH, token](const VqeConvergencePoint &pt) {
          if (token->cancelled.load()) {
            throw std::runtime_error("VQE optimization cancelled by user");
          }
          UiMessage msg;
          msg.type = UiMessage::Type::Convergence;
          msg.id = token->id;
          msg.convergence = pt;
          msg.molecule = isHeH ? QStringLiteral("HeH+") : QStringLiteral("H2");
          msg.text = QStringLiteral("%1 VQE Eval %2: E = %3 Ha%4")
                         .arg(msg.molecule)
                         .arg(pt.evaluation)
                         .arg(pt.totalEnergy, 0, 'f', 6)
                         .arg(pt.isBest ? QStringLiteral(" (New Best)") : QStringLiteral(""));
          post(msg);
        };

        if (!isHeH) {
          const auto res = optimizeH2(callback);
          UiMessage msg;
          msg.type = UiMessage::Type::Complete;
          msg.id = token->id;
          msg.text = QStringLiteral("H2 VQE Converged: Total E = %1 Ha (%2 evaluations)")
                         .arg(res.totalEnergy, 0, 'f', 6)
                         .arg(res.energyEvaluations);
          post(msg);
        } else {
          const auto res = optimizeHeH(callback);
          UiMessage msg;
          msg.type = UiMessage::Type::Complete;
          msg.id = token->id;
          msg.text = QStringLiteral("HeH+ VQE Converged: Total E = %1 Ha (%2 evaluations)")
                         .arg(res.totalEnergy, 0, 'f', 6)
                         .arg(res.energyEvaluations);
          post(msg);
        }
      } catch (const std::exception &err) {
        UiMessage msg;
        msg.type = token->cancelled.load() ? UiMessage::Type::Cancelled : UiMessage::Type::Failed;
        msg.id = token->id;
        msg.text = QString::fromStdString(err.what());
        post(msg);
      }
    });
  }

  void loadPresetExample(int index) {
    QString code;
    switch (index) {
    case 0:
      code = QStringLiteral(
          "# Quantum Bell State Preparation\n"
          "# Prepares maximally entangled state (|00⟩ + |11⟩) / √2\n"
          "qubit q0, q1\n"
          "H q0\n"
          "CNOT q0, q1\n"
          "measure q0, q1\n");
      break;
    case 1:
      code = QStringLiteral(
          "# Single-Qubit Superposition & Phase Rotation\n"
          "qubit q0\n"
          "H q0\n"
          "let theta = 3.141592653589793 / 4.0\n"
          "RZ(theta) q0\n"
          "H q0\n"
          "measure q0\n");
      break;
    case 2:
      code = QStringLiteral(
          "# 3-Qubit Greenberger-Horne-Zeilinger (GHZ) State\n"
          "# Prepares (|000⟩ + |111⟩) / √2\n"
          "qubit q0, q1, q2\n"
          "H q0\n"
          "CNOT q0, q1\n"
          "CNOT q1, q2\n"
          "measure q0, q1, q2\n");
      break;
    case 3:
      code = QStringLiteral(
          "# Deterministic Pauli-X Bit-Flip\n"
          "qubit q0\n"
          "X q0\n"
          "measure q0\n");
      break;
    case 4:
      code = formatHeHDsl(0.068498646886667, 4.976298509698156);
      break;
    case 5:
      code = QStringLiteral(
          "# Compiler Diagnostic Reporting Demonstration\n"
          "qubit q0\n"
          "H missing_qubit\n"
          "measure q0, q0\n"
          "let bad = 1.0 / 0.0\n");
      break;
    default:
      return;
    }
    isCompiling_ = true;
    if (dslEditor_) {
      dslEditor_->setPlainText(code);
    }
    isCompiling_ = false;
    compileDsl();
  }

  void navigateToDiagnostic(const compiler::Diagnostic &diagnostic) {
    if (!dslEditor_ || !diagnostic.sourceSpan.isValid()) return;

    const std::string source = dslEditor_->toPlainText().toStdString();
    const SourceSpan &span = diagnostic.sourceSpan;
    if (span.start.offset > source.size() ||
        span.end.offset > source.size() ||
        span.start.offset > static_cast<std::size_t>(
                                std::numeric_limits<qsizetype>::max()) ||
        span.end.offset > static_cast<std::size_t>(
                              std::numeric_limits<qsizetype>::max())) {
      return;
    }

    const auto documentPosition = [&source](std::size_t byteOffset) {
      return QString::fromUtf8(
                 source.data(), static_cast<qsizetype>(byteOffset))
          .size();
    };
    const qsizetype start = documentPosition(span.start.offset);
    const qsizetype end = documentPosition(span.end.offset);
    const int maximumPosition =
        std::max(0, dslEditor_->document()->characterCount() - 1);
    if (start < 0 || end < start || end > maximumPosition ||
        end > std::numeric_limits<int>::max()) {
      return;
    }

    QTextCursor cursor(dslEditor_->document());
    cursor.setPosition(static_cast<int>(start));
    cursor.setPosition(static_cast<int>(end), QTextCursor::KeepAnchor);
    dslEditor_->setTextCursor(cursor);
    dslEditor_->setFocus();
  }

  void navigateToDiagnosticRow(int row) {
    if (dslDebounceTimer_.isActive() || row < 0 ||
        row >= static_cast<int>(lastCompilation_.diagnostics.size())) {
      return;
    }
    navigateToDiagnostic(
        lastCompilation_.diagnostics[static_cast<std::size_t>(row)]);
  }

  void compileDsl() {
    if (!dslEditor_) return;
    const std::string source = dslEditor_->toPlainText().toStdString();
    compiler::CompilerOptions options;
    options.optimize = optimizeCheck_ ? optimizeCheck_->isChecked() : true;

    lastCompilation_ = compiler::compileProgram(source, options);

    if (lastCompilation_.success) {
      compileStatusLabel_->setText(
          QStringLiteral("COMPILED · %1 gate(s) · %2 logical layer(s) · %3 qubit(s)")
              .arg(lastCompilation_.instructions.size())
              .arg(lastCompilation_.optimizedDepth)
              .arg(lastCompilation_.qubitCount));
      compileStatusLabel_->setStyleSheet(QStringLiteral(
          "font-family:monospace;font-size:11px;font-weight:700;color:#38bdf8;padding:4px 8px;background:#0f172a;border:1px solid #1e293b;border-radius:2px"));
      diagnosticsTable_->setRowCount(0);
      diagnosticsTable_->hide();
      refreshCircuitFromCompilation();
    } else {
      compileStatusLabel_->setText(
          QStringLiteral("COMPILATION FAILED · %1 diagnostic(s)")
              .arg(lastCompilation_.diagnostics.size()));
      compileStatusLabel_->setStyleSheet(QStringLiteral(
          "font-family:monospace;font-size:11px;font-weight:700;color:#f43f5e;padding:4px 8px;background:#1e1014;border:1px solid #4c0519;border-radius:2px"));

      diagnosticsTable_->setRowCount(
          static_cast<int>(lastCompilation_.diagnostics.size()));
      for (std::size_t i = 0; i < lastCompilation_.diagnostics.size(); ++i) {
        const auto &diag = lastCompilation_.diagnostics[i];
        auto *lineItem = new QTableWidgetItem(
            diag.line > 0 ? QString::number(diag.line) : QStringLiteral("—"));
        auto *columnItem = new QTableWidgetItem(
            diag.column > 0 ? QString::number(diag.column)
                            : QStringLiteral("—"));
        lineItem->setTextAlignment(Qt::AlignCenter);
        columnItem->setTextAlignment(Qt::AlignCenter);
        auto *sevItem = new QTableWidgetItem(QStringLiteral("Error"));
        sevItem->setTextAlignment(Qt::AlignCenter);
        sevItem->setForeground(QColor(244, 63, 94));
        auto *msgItem = new QTableWidgetItem(QString::fromStdString(diag.message));
        msgItem->setForeground(QColor(254, 205, 211));

        diagnosticsTable_->setItem(static_cast<int>(i), 0, lineItem);
        diagnosticsTable_->setItem(static_cast<int>(i), 1, columnItem);
        diagnosticsTable_->setItem(static_cast<int>(i), 2, sevItem);
        diagnosticsTable_->setItem(static_cast<int>(i), 3, msgItem);
      }
      diagnosticsTable_->show();

      clearCircuitSimulationPresentation(
          QStringLiteral("Compilation error: simulation disabled"),
          QStringLiteral("Register: 0 Qubits"),
          QStringLiteral("STEP 0 / 0"));
      if (circuitView_) circuitView_->setCircuit({}, 0);
    }
  }

  void refreshCircuitFromCompilation() {
    if (!lastCompilation_.success) return;
    const bool optimizedView = viewOptimizedRadio_ && viewOptimizedRadio_->isChecked();
    const auto &instructions = optimizedView
        ? lastCompilation_.instructions
        : lastCompilation_.rawInstructions;
    const std::size_t displayedDepth = optimizedView
        ? lastCompilation_.optimizedDepth
        : lastCompilation_.rawDepth;
    if (compileStatusLabel_) {
      compileStatusLabel_->setText(
          QStringLiteral("COMPILED · %1 %2 gate(s) · %3 logical layer(s) · %4 qubit(s)")
              .arg(optimizedView ? QStringLiteral("optimized") : QStringLiteral("raw"))
              .arg(instructions.size())
              .arg(displayedDepth)
              .arg(lastCompilation_.qubitCount));
    }
    const int qubitCount = static_cast<int>(lastCompilation_.qubitCount);
    const QString metrics = QStringLiteral("Register: %1 Qubits · Simulation unavailable")
                                .arg(qubitCount);
    const QString stepText = QStringLiteral("STEP 0 / %1").arg(instructions.size());

    circuitPlayTimer_.stop();
    circuitSession_.reset();
    circuitAvailableStepLimit_ = 0;
    if (circuitView_) {
      circuitView_->setCircuit(instructions, qubitCount);
      circuitView_->setMinimumSize(circuitView_->sizeHint());
      circuitView_->setCurrentStep(0);
      if (circuitScrollArea_) {
        const int contentHeight = circuitView_->sizeHint().height();
        circuitScrollArea_->setMinimumHeight(
            std::clamp(contentHeight + 2, 180, 340));
      }
    }
    clearCircuitSimulationPresentation(
        QStringLiteral("Initializing simulation…"),
        QStringLiteral("Register: %1 Qubits · Simulation initializing").arg(qubitCount),
        stepText);

    try {
      if (smokeTest_ && injectNextSessionInitializationFailureForSmoke_) {
        injectNextSessionInitializationFailureForSmoke_ = false;
        throw std::length_error("Test-injected simulator initialization resource limit");
      }
      circuitSession_ = std::make_unique<SimulationSession>(instructions, qubitCount);
    } catch (const std::exception &error) {
      clearCircuitSimulationPresentation(
          QStringLiteral("Simulation unavailable for compiled circuit"),
          metrics,
          stepText);
      setCircuitFailureStatus(QStringLiteral("SIMULATOR"), error);
      return;
    }

    circuitAvailableStepLimit_ = circuitSession_->totalSteps();
    updateCircuitExecutionDisplay();
  }

  void clearCircuitSimulationPresentation(
      const QString &stageText,
      const QString &metricsText,
      const QString &stepText
  ) {
    circuitSession_.reset();
    circuitAvailableStepLimit_ = 0;
    circuitPlayTimer_.stop();
    if (circuitStepSlider_) {
      circuitStepSlider_->blockSignals(true);
      circuitStepSlider_->setRange(0, 0);
      circuitStepSlider_->setValue(0);
      circuitStepSlider_->blockSignals(false);
      circuitStepSlider_->setEnabled(false);
    }
    if (circuitStepLabel_) circuitStepLabel_->setText(stepText);
    if (circuitStageLabel_) circuitStageLabel_->setText(stageText);
    if (circuitMetricsLabel_) circuitMetricsLabel_->setText(metricsText);
    if (circuitMeasurementBanner_) {
      circuitMeasurementBanner_->clear();
      circuitMeasurementBanner_->hide();
    }
    if (circuitStateTable_) circuitStateTable_->setRowCount(0);
    if (circuitViewportHint_) circuitViewportHint_->hide();
    if (circuitResetBtn_) circuitResetBtn_->setEnabled(false);
    if (circuitStepBackBtn_) circuitStepBackBtn_->setEnabled(false);
    if (circuitStepFwdBtn_) circuitStepFwdBtn_->setEnabled(false);
    if (circuitPlayBtn_) circuitPlayBtn_->setEnabled(false);
    if (circuitPauseBtn_) circuitPauseBtn_->setEnabled(false);
  }

  QString circuitFailureMessage(
      const QString &operation,
      const std::exception &error
  ) const {
    const QString detail = QString::fromLocal8Bit(error.what());
    if (dynamic_cast<const std::length_error *>(&error) != nullptr ||
        dynamic_cast<const std::bad_alloc *>(&error) != nullptr) {
      return operation + QStringLiteral(" RESOURCE LIMIT: ") + detail;
    }
    if (dynamic_cast<const std::invalid_argument *>(&error) != nullptr ||
        dynamic_cast<const std::out_of_range *>(&error) != nullptr) {
      return operation + QStringLiteral(" INPUT ERROR: ") + detail;
    }
    return operation + QStringLiteral(" FAILED: ") + detail;
  }

  void setCircuitFailureStatus(
      const QString &operation,
      const std::exception &error
  ) {
    if (!compileStatusLabel_) return;
    compileStatusLabel_->setText(circuitFailureMessage(operation, error));
    compileStatusLabel_->setStyleSheet(QStringLiteral("color:#f43f5e;font-weight:bold;"));
  }

  void followCircuitOperation(std::size_t step) {
    if (!circuitView_ || !circuitScrollArea_ ||
        step >= circuitView_->instructionCount()) {
      if (circuitViewportHint_) circuitViewportHint_->hide();
      return;
    }
    QWidget *viewport = circuitScrollArea_->viewport();
    const QRect bounds = circuitView_->operationBounds(step);
    const QSize viewportSize = viewport->size();
    if (bounds.isEmpty() || viewportSize.isEmpty()) return;

    const auto targetOffset = [](int relativeBegin, int length,
                                 int viewportLength, int current, int maximum) {
      int target = current;
      if (length > viewportLength) {
        target = current + relativeBegin + length / 2 - viewportLength / 2;
      } else if (relativeBegin < 0) {
        target = current + relativeBegin;
      } else if (relativeBegin + length > viewportLength) {
        target = current + relativeBegin + length - viewportLength;
      }
      return std::clamp(target, 0, maximum);
    };
    const QPoint position = circuitView_->mapTo(viewport, bounds.topLeft());
    const int horizontal = targetOffset(
        position.x(), bounds.width(), viewportSize.width(),
        circuitScrollArea_->horizontalScrollBar()->value(),
        circuitScrollArea_->horizontalScrollBar()->maximum());
    const int vertical = targetOffset(
        position.y(), bounds.height(), viewportSize.height(),
        circuitScrollArea_->verticalScrollBar()->value(),
        circuitScrollArea_->verticalScrollBar()->maximum());
    circuitScrollArea_->horizontalScrollBar()->setValue(horizontal);
    circuitScrollArea_->verticalScrollBar()->setValue(vertical);

    const bool oversized = bounds.width() > viewportSize.width() ||
                           bounds.height() > viewportSize.height();
    if (circuitViewportHint_) {
      if (oversized) {
        circuitViewportHint_->setText(QStringLiteral(
            "Active operation exceeds the viewport; it is centered. Scroll to inspect its full control/target span."));
        circuitViewportHint_->show();
      } else {
        circuitViewportHint_->hide();
      }
    }
  }

  bool eventFilter(QObject *watched, QEvent *event) override {
    if (circuitScrollArea_ && watched == circuitScrollArea_->viewport() &&
        event->type() == QEvent::Resize) {
      QTimer::singleShot(0, this, [this] {
        if (circuitSession_) followCircuitOperation(circuitSession_->currentStep());
      });
    }
    return QMainWindow::eventFilter(watched, event);
  }

  void updateCircuitExecutionDisplay() {
    if (!circuitSession_) return;
    const std::size_t step = circuitSession_->currentStep();
    if (circuitStepSlider_) {
      circuitStepSlider_->blockSignals(true);
      circuitStepSlider_->setRange(0, static_cast<int>(circuitAvailableStepLimit_));
      circuitStepSlider_->setValue(static_cast<int>(step));
      circuitStepSlider_->blockSignals(false);
      circuitStepSlider_->setEnabled(circuitAvailableStepLimit_ > 0);
    }
    if (circuitView_) {
      circuitView_->setCurrentStep(step);
      followCircuitOperation(step);
    }
    refreshCircuitStateInspection();
  }

  void handleCircuitExecutionFailure(const std::exception &error) {
    circuitPlayTimer_.stop();
    if (circuitSession_) {
      circuitAvailableStepLimit_ = circuitSession_->currentStep();
      updateCircuitExecutionDisplay();
    }
    setCircuitFailureStatus(QStringLiteral("SIMULATION"), error);
  }

  void circuitPlay() {
    if (!circuitSession_) return;
    try {
      if (circuitSession_->isComplete() &&
          circuitAvailableStepLimit_ == circuitSession_->totalSteps()) {
        circuitSession_->reset();
        updateCircuitExecutionDisplay();
      }
      if (circuitSession_->currentStep() >= circuitAvailableStepLimit_) {
        circuitPlayTimer_.stop();
        refreshCircuitStateInspection();
        return;
      }
      circuitPlayTimer_.start();
      refreshCircuitStateInspection();
    } catch (const std::exception &error) {
      handleCircuitExecutionFailure(error);
    }
  }

  void circuitPause() {
    circuitPlayTimer_.stop();
    if (circuitSession_) refreshCircuitStateInspection();
  }

  void circuitStepForward() {
    if (!circuitSession_) return;
    if (circuitSession_->currentStep() >= circuitAvailableStepLimit_) {
      circuitPlayTimer_.stop();
      refreshCircuitStateInspection();
      return;
    }
    try {
      if (smokeTest_ && injectNextCircuitStepResourceFailureForSmoke_) {
        injectNextCircuitStepResourceFailureForSmoke_ = false;
        throw std::length_error("Test-injected snapshot memory limit");
      }
      circuitSession_->stepForward();
    } catch (const std::exception &error) {
      handleCircuitExecutionFailure(error);
      return;
    }
    if (circuitSession_->isComplete()) circuitPlayTimer_.stop();
    updateCircuitExecutionDisplay();
  }

  void circuitStepBackward() {
    if (!circuitSession_ || circuitSession_->currentStep() == 0) return;
    circuitPlayTimer_.stop();
    try {
      circuitSession_->stepBackward();
    } catch (const std::exception &error) {
      handleCircuitExecutionFailure(error);
      return;
    }
    updateCircuitExecutionDisplay();
  }

  void circuitReset() {
    if (!circuitSession_) return;
    circuitPlayTimer_.stop();
    try {
      circuitSession_->reset();
    } catch (const std::exception &error) {
      handleCircuitExecutionFailure(error);
      return;
    }
    updateCircuitExecutionDisplay();
  }

  void circuitSeek(int step) {
    if (!circuitSession_ || step < 0) return;
    circuitPlayTimer_.stop();
    const auto boundedStep = std::min(
        static_cast<std::size_t>(step), circuitAvailableStepLimit_);
    try {
      circuitSession_->seek(boundedStep);
    } catch (const std::exception &error) {
      handleCircuitExecutionFailure(error);
      return;
    }
    updateCircuitExecutionDisplay();
  }

  void refreshCircuitStateInspection() {
    if (!circuitSession_) return;

    const std::size_t step = circuitSession_->currentStep();
    const std::size_t total = circuitSession_->totalSteps();
    if (circuitStepLabel_) {
      circuitStepLabel_->setText(QStringLiteral("STEP %1 / %2").arg(step).arg(total));
    }

    if (circuitStageLabel_) {
      if (step == 0) {
        circuitStageLabel_->setText(QStringLiteral("Selected Stage: Initial |0...0⟩ state before instruction 1"));
      } else {
        const auto &instr = circuitSession_->instructions()[step - 1];
        circuitStageLabel_->setText(QStringLiteral("Selected Stage: After instruction %1 of %2 (%3)")
                                        .arg(step)
                                        .arg(total)
                                        .arg(QString::fromStdString(instr.gateName)));
      }
    }

    const auto measOpt = circuitSession_->lastMeasurementResult();
    if (circuitMeasurementBanner_) {
      if (measOpt.has_value()) {
        circuitMeasurementBanner_->setText(QStringLiteral(
            "PROJECTIVE MEASUREMENT COLLAPSED: Measured outcome is '%1' · State projected and renormalized")
                                               .arg(QString::fromStdString(*measOpt)));
        circuitMeasurementBanner_->show();
      } else {
        circuitMeasurementBanner_->hide();
      }
    }

    const auto &state = circuitSession_->currentState();
    const int n = circuitSession_->qubitCount();
    const std::size_t dimension = state.size();

    if (circuitStateTable_) {
      circuitStateTable_->setRowCount(0);
      constexpr std::size_t kMaxVisibleRows = 128;
      const bool truncate = (dimension > kMaxVisibleRows);
      const std::size_t rowsToDisplay = truncate ? kMaxVisibleRows : dimension;

      circuitStateTable_->setRowCount(static_cast<int>(rowsToDisplay + (truncate ? 1 : 0)));

      for (std::size_t i = 0; i < rowsToDisplay; ++i) {
        const auto &amp = state[i];
        const double p = std::norm(amp);
        const QString basisStr = QStringLiteral("|%1⟩ (index %2)")
                                     .arg(QString::number(i, 2).rightJustified(n, QLatin1Char('0')))
                                     .arg(i);
        auto *basisItem = new QTableWidgetItem(basisStr);
        basisItem->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);

        const QString ampStr = QStringLiteral("%1 %2 %3 i")
                                   .arg(amp.real(), 0, 'f', 6)
                                   .arg(amp.imag() < 0.0 ? QStringLiteral("−") : QStringLiteral("+"))
                                   .arg(std::abs(amp.imag()), 0, 'f', 6);
        auto *ampItem = new QTableWidgetItem(ampStr);
        ampItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);

        auto *probItem = new QTableWidgetItem(QString::number(p, 'f', 6));
        probItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);

        auto *barWidget = new QWidget;
        auto *barLayout = new QHBoxLayout(barWidget);
        barLayout->setContentsMargins(4, 2, 4, 2);
        auto *prog = new QProgressBar;
        prog->setRange(0, 1000);
        prog->setValue(static_cast<int>(p * 1000.0));
        prog->setTextVisible(false);
        prog->setFixedHeight(10);
        barLayout->addWidget(prog);

        circuitStateTable_->setItem(static_cast<int>(i), 0, basisItem);
        circuitStateTable_->setItem(static_cast<int>(i), 1, ampItem);
        circuitStateTable_->setItem(static_cast<int>(i), 2, probItem);
        circuitStateTable_->setCellWidget(static_cast<int>(i), 3, barWidget);
      }

      if (truncate) {
        auto *moreItem = new QTableWidgetItem(
            QStringLiteral("... and %1 additional computational basis states").arg(dimension - kMaxVisibleRows));
        moreItem->setForeground(QColor(100, 116, 139));
        circuitStateTable_->setItem(static_cast<int>(kMaxVisibleRows), 0, moreItem);
      }
    }

    const double norm = circuitSession_->stateNorm();
    if (circuitMetricsLabel_) {
      circuitMetricsLabel_->setText(
          QStringLiteral("Register: %1 Qubits · Dimension: %2 · Norm: %3 (%4)")
              .arg(n)
              .arg(dimension)
              .arg(norm, 0, 'f', 12)
              .arg(std::abs(norm - 1.0) < 1e-10 ? QStringLiteral("Normalized") : QStringLiteral("Invalid")));
    }

    if (circuitStepSlider_) circuitStepSlider_->setEnabled(circuitAvailableStepLimit_ > 0);
    if (circuitStepBackBtn_) circuitStepBackBtn_->setEnabled(step > 0);
    if (circuitStepFwdBtn_) circuitStepFwdBtn_->setEnabled(step < circuitAvailableStepLimit_);
    if (circuitResetBtn_) circuitResetBtn_->setEnabled(step > 0);
    if (circuitPlayBtn_) {
      const bool canRestart = step == total && circuitAvailableStepLimit_ == total && total > 0;
      circuitPlayBtn_->setEnabled(step < circuitAvailableStepLimit_ || canRestart);
    }
    if (circuitPauseBtn_) circuitPauseBtn_->setEnabled(circuitPlayTimer_.isActive());
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
          UiMessage msgProg;
          msgProg.type = UiMessage::Type::Progress;
          msgProg.id = token->id;
          msgProg.completed = index;
          msgProg.total = total;
          msgProg.sweep = sweep;
          msgProg.text = QStringLiteral("Preparing sample %1/%2 at %3 Å")
                             .arg(index + 1).arg(total).arg(distance, 0, 'f', 3);
          post(msgProg);
          const ProgressCallback callback{[this, token, index, total, sweep](const QString &text) {
            UiMessage cbMsg;
            cbMsg.type = UiMessage::Type::Progress;
            cbMsg.id = token->id;
            cbMsg.completed = index;
            cbMsg.total = total;
            cbMsg.sweep = sweep;
            cbMsg.text = text;
            post(cbMsg);
          }};
          EnergyResult result = calculateAtGeometry(root, distance, token, callback);
          if (token->cancelled.load()) throw std::runtime_error("Calculation cancelled");
          UiMessage msgPoint;
          msgPoint.type = UiMessage::Type::Point;
          msgPoint.id = token->id;
          msgPoint.completed = index + 1;
          msgPoint.total = total;
          msgPoint.sweep = sweep;
          msgPoint.text = QStringLiteral("Completed actual calculation at %1 Å")
                              .arg(result.distanceAngstrom, 0, 'f', 3);
          msgPoint.result = result;
          post(msgPoint);
        }
        UiMessage msgDone;
        msgDone.type = UiMessage::Type::Complete;
        msgDone.id = token->id;
        msgDone.completed = total;
        msgDone.total = total;
        msgDone.sweep = sweep;
        msgDone.text = QStringLiteral("All requested calculations completed.");
        post(msgDone);
      } catch (const std::exception &error) {
        const bool cancelled = token->cancelled.load() ||
            QString::fromUtf8(error.what()).contains(QStringLiteral("cancelled"), Qt::CaseInsensitive);
        UiMessage msgErr;
        msgErr.type = cancelled ? UiMessage::Type::Cancelled : UiMessage::Type::Failed;
        msgErr.id = token->id;
        msgErr.sweep = sweep;
        msgErr.text = QString::fromUtf8(error.what());
        post(msgErr);
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
      case UiMessage::Type::Convergence:
        if (convergencePlot_) {
          convergencePlot_->addPoint(message.convergence);
        }
        if (vqeStatusLabel_) {
          vqeStatusLabel_->setText(message.text);
        }
        break;
      case UiMessage::Type::Point:
        plotResults_.push_back(message.result);
        if (message.sweep) currentSweepResults_.push_back(message.result);
        updateEnergyPlotSeries();
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
        if (liveH2VqeBtn_) liveH2VqeBtn_->setEnabled(true);
        if (liveHeHVqeBtn_) liveHeHVqeBtn_->setEnabled(true);
        if (cancelVqeBtn_) cancelVqeBtn_->setEnabled(false);
        progress_->setValue(100);
        progressLabel_->setText(message.text);
        if (vqeStatusLabel_) vqeStatusLabel_->setText(message.text);
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
        if (liveH2VqeBtn_) liveH2VqeBtn_->setEnabled(true);
        if (liveHeHVqeBtn_) liveHeHVqeBtn_->setEnabled(true);
        if (cancelVqeBtn_) cancelVqeBtn_->setEnabled(false);
        progressLabel_->setText((message.type == UiMessage::Type::Failed ? QStringLiteral("Failed: ") : QStringLiteral("Cancelled: ")) + message.text);
        if (vqeStatusLabel_) vqeStatusLabel_->setText((message.type == UiMessage::Type::Failed ? QStringLiteral("VQE Error: ") : QStringLiteral("VQE Cancelled: ")) + message.text);
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

    // 1. Verify Circuit Workbench interactive compilation and stepping
    if (!dslEditor_ || !circuitSession_) {
      finishSmokeTest(false, QStringLiteral("circuit workbench not initialized"));
      return;
    }
    if (circuitSession_->qubitCount() != 2 || circuitSession_->totalSteps() != 3) {
      finishSmokeTest(false, QStringLiteral("Bell state preset did not compile to 2 qubits, 3 steps"));
      return;
    }
    // Verify initial state |00> = 1.0
    if (std::abs(circuitSession_->currentState()[0] - std::complex<double>(1.0, 0.0)) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("initial state is not |00>"));
      return;
    }
    // Step forward 1: H q0
    circuitStepFwdBtn_->click();
    if (circuitSession_->currentStep() != 1) {
      finishSmokeTest(false, QStringLiteral("step forward 1 did not advance step to 1"));
      return;
    }
    const double invSqrt2 = 1.0 / std::sqrt(2.0);
    if (std::abs(circuitSession_->currentState()[0] - std::complex<double>(invSqrt2, 0.0)) > 1e-12 ||
        std::abs(circuitSession_->currentState()[1] - std::complex<double>(invSqrt2, 0.0)) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("step 1 state does not match Hadamard superposition"));
      return;
    }
    // Step forward 2: CNOT q0, q1
    circuitStepFwdBtn_->click();
    if (circuitSession_->currentStep() != 2) {
      finishSmokeTest(false, QStringLiteral("step forward 2 did not advance step to 2"));
      return;
    }
    const auto bellProbs = circuitSession_->basisProbabilities();
    if (std::abs(bellProbs[0] - 0.5) > 1e-12 || std::abs(bellProbs[3] - 0.5) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("step 2 state does not match Bell state (|00> + |11>)/sqrt(2)"));
      return;
    }
    // Step forward 3: Measure q0, q1
    circuitStepFwdBtn_->click();
    if (!circuitSession_->isComplete() || !circuitSession_->lastMeasurementResult().has_value()) {
      finishSmokeTest(false, QStringLiteral("step 3 measurement did not collapse or record outcome"));
      return;
    }
    const std::string bellMeas = *circuitSession_->lastMeasurementResult();
    if (bellMeas != "00" && bellMeas != "11") {
      finishSmokeTest(false, QStringLiteral("Bell state measurement outcome was not '00' or '11'"));
      return;
    }
    // Step backward: verify superposition restored
    circuitStepBackBtn_->click();
    if (circuitSession_->currentStep() != 2 ||
        std::abs(circuitSession_->basisProbabilities()[0] - 0.5) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("step backward did not restore unmeasured Bell superposition"));
      return;
    }
    // Step forward: verify measurement outcome preserved
    circuitStepFwdBtn_->click();
    if (*circuitSession_->lastMeasurementResult() != bellMeas) {
      finishSmokeTest(false, QStringLiteral("step forward did not preserve recorded measurement outcome"));
      return;
    }
    if (!circuitMeasurementBanner_ || circuitMeasurementBanner_->isHidden()) {
      finishSmokeTest(false, QStringLiteral("completed Bell run did not expose its measurement result"));
      return;
    }

    auto compileSmokeSource = [this](const QString &source) {
      isCompiling_ = true;
      dslEditor_->setPlainText(source);
      isCompiling_ = false;
      compileDsl();
    };

    QStringList unsupportedQubits;
    for (int i = 0; i < 25; ++i) {
      unsupportedQubits << QStringLiteral("q%1").arg(i);
    }
    compileSmokeSource(QStringLiteral("qubit %1\nH q0\nmeasure q0\n")
                           .arg(unsupportedQubits.join(QStringLiteral(", "))));
    if (!lastCompilation_.success || circuitSession_ ||
        circuitView_->qubitCount() != 25 || circuitView_->instructionCount() != 2 ||
        circuitView_->currentStep() != 0 ||
        !compileStatusLabel_->text().contains(QStringLiteral("SIMULATOR INPUT ERROR")) ||
        !compileStatusLabel_->text().contains(QStringLiteral("24 qubits")) ||
        !circuitMeasurementBanner_->isHidden() || circuitStateTable_->rowCount() != 0 ||
        circuitStepSlider_->isEnabled() || circuitStepFwdBtn_->isEnabled() ||
        circuitPlayBtn_->isEnabled()) {
      finishSmokeTest(false, QStringLiteral("failed simulator initialization retained stale output or enabled execution"));
      return;
    }

    loadPresetExample(0);
    circuitStepFwdBtn_->click();
    circuitStepFwdBtn_->click();
    circuitStepFwdBtn_->click();
    if (!circuitSession_->lastMeasurementResult().has_value() ||
        circuitMeasurementBanner_->isHidden()) {
      finishSmokeTest(false, QStringLiteral("Bell measurement setup for compilation-failure check did not complete"));
      return;
    }
    compileSmokeSource(QStringLiteral(
        "qubit q0\nmeasure q0, q0\nH missing\n"));
    if (lastCompilation_.success || circuitSession_ ||
        circuitView_->qubitCount() != 0 || circuitView_->instructionCount() != 0 ||
        !circuitMeasurementBanner_->isHidden() || circuitStateTable_->rowCount() != 0 ||
        diagnosticsTable_->isHidden() ||
        !compileStatusLabel_->text().contains(QStringLiteral("COMPILATION FAILED"))) {
      finishSmokeTest(false, QStringLiteral("compilation failure did not preserve diagnostics while clearing old run output"));
      return;
    }

    const auto clickDiagnosticRow = [this](int row) {
      QWidget *viewport = diagnosticsTable_->viewport();
      const QRect cell = diagnosticsTable_->visualRect(
          diagnosticsTable_->model()->index(row, 0));
      if (!cell.isValid()) return false;
      const QPointF localPosition(cell.center());
      const QPointF globalPosition(viewport->mapToGlobal(cell.center()));
      QMouseEvent press(QEvent::MouseButtonPress, localPosition, globalPosition,
                        Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
      QApplication::sendEvent(viewport, &press);
      QMouseEvent release(QEvent::MouseButtonRelease, localPosition, globalPosition,
                          Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
      QApplication::sendEvent(viewport, &release);
      return true;
    };
    if (lastCompilation_.diagnostics.size() != 2 ||
        diagnosticsTable_->columnCount() != 4 ||
        !diagnosticsTable_->item(0, 0) ||
        !diagnosticsTable_->item(0, 1) ||
        diagnosticsTable_->item(0, 0)->text() != QStringLiteral("2") ||
        diagnosticsTable_->item(0, 1)->text() != QStringLiteral("13") ||
        !clickDiagnosticRow(0)) {
      finishSmokeTest(false, QStringLiteral("diagnostic table did not expose precise source coordinates"));
      return;
    }
    const QTextCursor duplicateCursor = dslEditor_->textCursor();
    if (duplicateCursor.selectedText() != QStringLiteral("q0") ||
        duplicateCursor.selectionStart() != 21 ||
        duplicateCursor.selectionEnd() != 23 ||
        lastCompilation_.diagnostics.size() != 2) {
      finishSmokeTest(false, QStringLiteral("clicking a diagnostic row did not select its exact duplicate token"));
      return;
    }
    if (!clickDiagnosticRow(1) ||
        dslEditor_->textCursor().selectedText() != QStringLiteral("missing") ||
        dslEditor_->textCursor().selectionStart() != 26) {
      finishSmokeTest(false, QStringLiteral("clicking a later diagnostic did not navigate to its source span"));
      return;
    }
    compileSmokeSource(QStringLiteral("#π\nH q0 @\n"));
    if (lastCompilation_.diagnostics.size() != 1 ||
        lastCompilation_.diagnostics[0].line != 2 ||
        lastCompilation_.diagnostics[0].column != 6 ||
        lastCompilation_.diagnostics[0].sourceSpan.start.offset != 9 ||
        lastCompilation_.diagnostics[0].sourceSpan.end.offset != 10 ||
        !diagnosticsTable_->item(0, 0) ||
        diagnosticsTable_->item(0, 0)->text() != QStringLiteral("2") ||
        diagnosticsTable_->item(0, 1)->text() != QStringLiteral("6") ||
        !clickDiagnosticRow(0) ||
        dslEditor_->textCursor().selectedText() != QStringLiteral("@") ||
        dslEditor_->textCursor().selectionStart() != 8 ||
        dslEditor_->textCursor().selectionEnd() != 9) {
      finishSmokeTest(false, QStringLiteral("UTF-8 byte offsets did not map to the correct Qt source character"));
      return;
    }
    const QString supplementaryCharacter =
        QString::fromUtf8("\xF0\x9F\x98\x80");
    compileSmokeSource(QStringLiteral("qubit q0\nH q0 ") +
                       supplementaryCharacter + QStringLiteral("\n"));
    if (lastCompilation_.diagnostics.size() != 1 ||
        lastCompilation_.diagnostics[0].line != 2 ||
        lastCompilation_.diagnostics[0].column != 6 ||
        lastCompilation_.diagnostics[0].sourceSpan.start.offset != 14 ||
        lastCompilation_.diagnostics[0].sourceSpan.end.offset != 18 ||
        lastCompilation_.diagnostics[0].sourceSpan.end.column != 10 ||
        !clickDiagnosticRow(0) ||
        dslEditor_->textCursor().selectedText() != supplementaryCharacter ||
        dslEditor_->textCursor().selectionStart() != 14 ||
        dslEditor_->textCursor().selectionEnd() != 16) {
      finishSmokeTest(false, QStringLiteral("supplementary UTF-8 diagnostic did not select the complete source character"));
      return;
    }


    compileSmokeSource(QStringLiteral("qubit q0\nH"));
    if (lastCompilation_.diagnostics.size() != 1 ||
        !lastCompilation_.diagnostics[0].sourceSpan.isValid() ||
        lastCompilation_.diagnostics[0].sourceSpan.start.offset != 10 ||
        lastCompilation_.diagnostics[0].sourceSpan.end.offset != 10 ||
        !clickDiagnosticRow(0) ||
        dslEditor_->textCursor().hasSelection() ||
        dslEditor_->textCursor().position() != 10) {
      finishSmokeTest(false, QStringLiteral("EOF diagnostic did not navigate to its zero-length insertion point"));
      return;
    }

    dslEditor_->setPlainText(QStringLiteral("qubit q0\nX q0\n"));
    const int cursorBeforeStaleClick = dslEditor_->textCursor().position();
    if (!dslDebounceTimer_.isActive() || !clickDiagnosticRow(0) ||
        dslEditor_->textCursor().position() != cursorBeforeStaleClick ||
        dslEditor_->textCursor().hasSelection()) {
      finishSmokeTest(false, QStringLiteral("a pending source edit allowed navigation using stale diagnostic offsets"));
      return;
    }
    dslDebounceTimer_.stop();
    compileSmokeSource(QStringLiteral("qubit q0\nH"));

    if (lastCompilation_.diagnostics.size() != 1 ||
        !lastCompilation_.diagnostics[0].sourceSpan.isValid() ||
        lastCompilation_.diagnostics[0].sourceSpan.start.offset != 10) {
      finishSmokeTest(false, QStringLiteral("EOF diagnostic setup failed after stale-location check"));
      return;
    }

    auto &unlocatedDiagnostic = lastCompilation_.diagnostics[0];
    unlocatedDiagnostic.line = 0;
    unlocatedDiagnostic.column = 0;
    unlocatedDiagnostic.sourceSpan = {};
    diagnosticsTable_->item(0, 0)->setText(QStringLiteral("—"));
    diagnosticsTable_->item(0, 1)->setText(QStringLiteral("—"));
    const int cursorBeforeUnlocatedClick = dslEditor_->textCursor().position();
    if (!clickDiagnosticRow(0) ||
        dslEditor_->textCursor().position() != cursorBeforeUnlocatedClick ||
        dslEditor_->textCursor().hasSelection()) {
      finishSmokeTest(false, QStringLiteral("unlocated diagnostic moved the editor to an arbitrary source position"));
      return;
    }

    compileSmokeSource(QStringLiteral("qubit q0\nH q0\n"));
    circuitStepFwdBtn_->click();
    if (!lastCompilation_.success || !circuitSession_ ||
        circuitStateTable_->rowCount() != 2 ||
        !circuitStateTable_->item(0, 0) || !circuitStateTable_->item(1, 0) ||
        !circuitStateTable_->item(0, 0)->text().contains(QStringLiteral("index 0")) ||
        !circuitStateTable_->item(1, 0)->text().contains(QStringLiteral("index 1")) ||
        !circuitStateTable_->item(0, 2) || !circuitStateTable_->item(1, 2) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 0.5) > 1e-12 ||
        std::abs(circuitStateTable_->item(1, 2)->text().toDouble() - 0.5) > 1e-12 ||
        !circuitMetricsLabel_->text().contains(QStringLiteral("Normalized"))) {
      finishSmokeTest(false, QStringLiteral("one-qubit probability display has incorrect range or values"));
      return;
    }

    QStringList largeQubits;
    for (int i = 0; i < 8; ++i) {
      largeQubits << QStringLiteral("q%1").arg(i);
    }
    compileSmokeSource(QStringLiteral("qubit %1\nH q0\n")
                           .arg(largeQubits.join(QStringLiteral(", "))));
    circuitStepFwdBtn_->click();
    if (!lastCompilation_.success || !circuitSession_ ||
        circuitStateTable_->rowCount() != 129 ||
        !circuitStateTable_->item(0, 0) || !circuitStateTable_->item(1, 0) ||
        !circuitStateTable_->item(0, 0)->text().contains(QStringLiteral("index 0")) ||
        !circuitStateTable_->item(1, 0)->text().contains(QStringLiteral("index 1")) ||
        !circuitStateTable_->item(0, 2) || !circuitStateTable_->item(1, 2) ||
        !circuitStateTable_->item(128, 0) ||
        !circuitStateTable_->item(128, 0)->text().contains(QStringLiteral("128 additional")) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 0.5) > 1e-12 ||
        std::abs(circuitStateTable_->item(1, 2)->text().toDouble() - 0.5) > 1e-12 ||
        !circuitMetricsLabel_->text().contains(QStringLiteral("Normalized"))) {
      finishSmokeTest(false, QStringLiteral("large-state inspector did not display a bounded 128-row prefix"));
      return;
    }

    loadPresetExample(0);
    if (!circuitSession_ || circuitSession_->currentStep() != 0) {
      finishSmokeTest(false, QStringLiteral("Bell preset did not restore a fresh interactive session"));
      return;
    }
    circuitStepFwdBtn_->click();
    if (circuitSession_->currentStep() != 1 || circuitStepSlider_->value() != 1 ||
        circuitView_->currentStep() != 1 || !circuitStateTable_->item(0, 2) ||
        !circuitStateTable_->item(1, 2) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 0.5) > 1e-12 ||
        std::abs(circuitStateTable_->item(1, 2)->text().toDouble() - 0.5) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("successful pre-failure step was not reflected consistently"));
      return;
    }
    injectNextCircuitStepResourceFailureForSmoke_ = true;
    circuitPlayBtn_->click();
    if (!circuitPlayTimer_.isActive()) {
      finishSmokeTest(false, QStringLiteral("playback did not start before injected step failure"));
      return;
    }
    QEventLoop failureLoop;
    QTimer::singleShot(circuitPlayTimer_.interval() + 50, &failureLoop, &QEventLoop::quit);
    failureLoop.exec();
    if (circuitPlayTimer_.isActive() || circuitSession_->currentStep() != 1 ||
        circuitStepSlider_->value() != 1 || circuitStepSlider_->maximum() != 1 ||
        circuitView_->currentStep() != 1 || circuitStepFwdBtn_->isEnabled() ||
        circuitPlayBtn_->isEnabled() ||
        !compileStatusLabel_->text().contains(QStringLiteral("SIMULATION RESOURCE LIMIT")) ||
        !compileStatusLabel_->text().contains(QStringLiteral("Test-injected snapshot memory limit")) ||
        !circuitStateTable_->item(0, 2) || !circuitStateTable_->item(1, 2) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 0.5) > 1e-12 ||
        std::abs(circuitStateTable_->item(1, 2)->text().toDouble() - 0.5) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("failed playback step advanced the UI or left its timer active"));
      return;
    }
    circuitPlay();
    circuitStepForward();
    if (circuitPlayTimer_.isActive() || circuitSession_->currentStep() != 1 ||
        circuitStepSlider_->value() != 1 || circuitView_->currentStep() != 1) {
      finishSmokeTest(false, QStringLiteral("repeated playback attempts changed the last valid step or restarted the timer"));
      return;
    }

    loadPresetExample(0);
    if (!circuitSession_ || circuitSession_->currentStep() != 0) {
      finishSmokeTest(false, QStringLiteral("Bell preset did not recover after resource failure"));
      return;
    }
    circuitStepFwdBtn_->click();
    circuitStepFwdBtn_->click();
    circuitStepSlider_->setValue(1);
    if (circuitSession_->currentStep() != 1 || circuitStepSlider_->value() != 1 ||
        circuitView_->currentStep() != 1 || !circuitStateTable_->item(0, 2) ||
        !circuitStateTable_->item(1, 2) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 0.5) > 1e-12 ||
        std::abs(circuitStateTable_->item(1, 2)->text().toDouble() - 0.5) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("slider seek did not synchronize the last selected Bell-state stage"));
      return;
    }
    circuitResetBtn_->click();
    if (circuitSession_->currentStep() != 0 || circuitStepSlider_->value() != 0 ||
        circuitView_->currentStep() != 0 || !circuitStateTable_->item(0, 2) ||
        std::abs(circuitStateTable_->item(0, 2)->text().toDouble() - 1.0) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("reset did not synchronize the initial circuit and state"));
      return;
    }

    // 1b. Save and reopen a stepped Bell run with attached saved chemistry metadata.
    const QString smokeQcomp = QStringLiteral("/tmp/smoke_exp.qcomp");
    const QString smokeSvg = QStringLiteral("/tmp/smoke_circ.svg");
    const QString smokePng = QStringLiteral("/tmp/smoke_circ.png");
    const QString smokeCsv = QStringLiteral("/tmp/smoke_conv.csv");
    const QString invalidOpenPath = QStringLiteral("/tmp/smoke_invalid_open.qcomp");
    const QString initFailurePath = QStringLiteral("/tmp/smoke_init_failure.qcomp");
    optimizeCheck_->setChecked(false);
    loadPresetExample(0);
    for (int step = 0; step < 3; ++step) circuitStepFwdBtn_->click();
    if (!circuitSession_ || circuitSession_->currentStep() != 3 ||
        !circuitSession_->lastMeasurementResult()) {
      finishSmokeTest(false, QStringLiteral("Bell persistence fixture did not reach its recorded measurement"));
      return;
    }

    EnergyResult savedResult;
    savedResult.distanceAngstrom = 0.735;
    savedResult.distanceBohr = 0.735 / 0.52917721092;
    savedResult.parameter = 0.123456789;
    savedResult.electronicEnergy = -1.857275030202;
    savedResult.exactElectronicEnergy = -1.857275030202;
    savedResult.referenceElectronicEnergy = -1.857275030202;
    savedResult.nuclearRepulsionEnergy = 0.719968994449;
    savedResult.totalEnergy = -1.137306035753;
    savedResult.exactTotalEnergy = -1.137306035753;
    savedResult.evaluations = 62;
    savedResult.converged = true;
    latestResult_ = savedResult;
    plotResults_ = {savedResult};
    currentSweepResults_ = {savedResult};
    updateResultFields(savedResult);
    updateEnergyPlotSeries();
    if (convergencePlot_) {
      convergencePlot_->reset(savedResult.exactTotalEnergy, QStringLiteral("H2"));
      convergencePlot_->addPoint({1, 0.0, 0.0, savedResult.electronicEnergy,
                                  savedResult.totalEnergy, true});
      convergencePlot_->addPoint({2, savedResult.parameter, 0.0, savedResult.electronicEnergy,
                                  savedResult.totalEnergy, true});
    }
    if (!saveCurrentExperiment(smokeQcomp) || !QFileInfo::exists(smokeQcomp)) {
      finishSmokeTest(false, QStringLiteral("saveCurrentExperiment failed to write .qcomp file"));
      return;
    }
    const ExperimentBundle savedBundle = loadExperimentBundle(smokeQcomp.toStdString());
    if (savedBundle.dslSource != dslEditor_->toPlainText().toStdString() ||
        savedBundle.optimize || savedBundle.currentStep != 3 ||
        savedBundle.recordedMeasurements.size() != 1 ||
        savedBundle.recordedMeasurements[0].second != *circuitSession_->lastMeasurementResult() ||
        !savedBundle.optimizerParameterRadians ||
        std::abs(*savedBundle.optimizerParameterRadians - savedResult.parameter) > 1e-15 ||
        !savedBundle.referenceElectronicEnergy ||
        savedBundle.convergenceTrajectory.size() != 2) {
      finishSmokeTest(false, QStringLiteral("saved .qcomp bundle omitted source, options, session, chemistry, or convergence data"));
      return;
    }

    // A rejected unsupported save must preserve the previously valid destination.
    compileSmokeSource(QStringLiteral("qubit %1\nH q0\nmeasure q0\n")
                           .arg(unsupportedQubits.join(QStringLiteral(", "))));
    const bool rejectedSave = !saveCurrentExperiment(smokeQcomp);
    const ExperimentBundle afterRejectedSave = loadExperimentBundle(smokeQcomp.toStdString());
    if (!rejectedSave ||
        !compileStatusLabel_->text().contains(QStringLiteral("SAVE FAILED")) ||
        afterRejectedSave.id != savedBundle.id ||
        afterRejectedSave.dslSource != savedBundle.dslSource ||
        afterRejectedSave.title != savedBundle.title) {
      finishSmokeTest(false, QStringLiteral("rejected save replaced or falsely reported success for the valid destination"));
      return;
    }

    if (!openExperimentFile(smokeQcomp) || !circuitSession_) {
      finishSmokeTest(false, QStringLiteral("openExperimentFile failed to restore the valid bundle"));
      return;
    }
    const auto &restoredState = circuitSession_->currentState();
    bool restoredStateMatches = restoredState.size() == savedBundle.statevector.size();
    for (std::size_t i = 0; restoredStateMatches && i < restoredState.size(); ++i) {
      restoredStateMatches = std::abs(restoredState[i] - savedBundle.statevector[i]) <= 1e-12;
    }
    bool inspectorMatches = circuitStateTable_->rowCount() == 4;
    for (std::size_t i = 0; inspectorMatches && i < savedBundle.basisProbabilities.size(); ++i) {
      const auto *probability = circuitStateTable_->item(static_cast<int>(i), 2);
      inspectorMatches = probability &&
          std::abs(probability->text().toDouble() - savedBundle.basisProbabilities[i]) <= 1e-6;
    }
    if (dslEditor_->toPlainText().toStdString() != savedBundle.dslSource ||
        optimizeCheck_->isChecked() || circuitSession_->currentStep() != 3 ||
        circuitStepSlider_->value() != 3 || circuitView_->currentStep() != 3 ||
        !circuitSession_->lastMeasurementResult() ||
        *circuitSession_->lastMeasurementResult() != savedBundle.recordedMeasurements[0].second ||
        !restoredStateMatches || !inspectorMatches ||
        !latestResult_ || std::abs(latestResult_->parameter - savedResult.parameter) > 1e-15 ||
        plotResults_.size() != 1 || !currentSweepResults_.empty() ||
        !convergencePlot_ || convergencePlot_->points().size() != 2 ||
        std::abs(convergencePlot_->exactReference() - savedResult.exactTotalEnergy) > 1e-15 ||
        !appStatus_->text().contains(QStringLiteral("not recalculated")) ||
        !circuitMeasurementBanner_->text().contains(QString::fromStdString(savedBundle.recordedMeasurements[0].second))) {
      finishSmokeTest(false, QStringLiteral("reopening .qcomp did not restore source, measurement, selected step, or saved-result metadata consistently"));
      return;
    }
    const auto &h2Series = energyPlot_->seriesList();
    if (h2Series.empty() || h2Series[0].points.size() != 1 ||
        std::abs(h2Series[0].points[0].totalEnergy - savedResult.totalEnergy) > 1e-15 ||
        !optimizerValue_->text().contains(QStringLiteral("θ=0.123456789 rad"))) {
      finishSmokeTest(false, QStringLiteral("saved chemistry summary was not restored as saved data"));
      return;
    }

    if (!exportCircuitToSvgFile(smokeSvg) || !QFileInfo::exists(smokeSvg)) {
      finishSmokeTest(false, QStringLiteral("exportCircuitToSvgFile failed to write SVG"));
      return;
    }
    if (!exportCircuitToPngFile(smokePng) || !QFileInfo::exists(smokePng)) {
      finishSmokeTest(false, QStringLiteral("exportCircuitToPngFile failed to write PNG"));
      return;
    }
    if (!exportConvergenceToCsvFile(smokeCsv) || !QFileInfo::exists(smokeCsv)) {
      finishSmokeTest(false, QStringLiteral("exportConvergenceToCsvFile failed to write CSV"));
      return;
    }
    QFile::remove(smokeSvg);
    QFile::remove(smokePng);
    QFile::remove(smokeCsv);

    // A source-compilation failure after parsing must not display OPENED or retain prior output.
    ExperimentBundle invalidSourceBundle;
    invalidSourceBundle.dslSource = "qubit q0\nH missing\n";
    saveExperimentBundle(invalidOpenPath.toStdString(), invalidSourceBundle);
    if (openExperimentFile(invalidOpenPath) || !compileStatusLabel_->text().contains(QStringLiteral("COMPILATION FAILED")) ||
        compileStatusLabel_->text().contains(QStringLiteral("OPENED")) ||
        dslEditor_->toPlainText() != QStringLiteral("qubit q0\nH missing\n") ||
        circuitSession_ || latestResult_ ||
        circuitStateTable_->rowCount() != 0 || !circuitMeasurementBanner_->isHidden() ||
        !convergencePlot_->points().empty()) {
      finishSmokeTest(false, QStringLiteral("failed compilation while opening a bundle left stale output or reported success"));
      return;
    }
    if (!openExperimentFile(smokeQcomp)) {
      finishSmokeTest(false, QStringLiteral("valid saved bundle did not reopen after a compilation failure"));
      return;
    }

    // Inject initialization failure at the actual compilation-to-session boundary.
    ExperimentBundle initFailureBundle;
    initFailureBundle.dslSource = "qubit q0\nH q0\nmeasure q0\n";
    initFailureBundle.qubitCount = 1;
    initFailureBundle.rawInstructions = {
        {InstructionType::ApplyGate, "H", {0}, 0.0},
        {InstructionType::Measure, "measure", {0}, 0.0}};
    initFailureBundle.instructions = initFailureBundle.rawInstructions;
    saveExperimentBundle(initFailurePath.toStdString(), initFailureBundle);
    injectNextSessionInitializationFailureForSmoke_ = true;
    if (openExperimentFile(initFailurePath) || !compileStatusLabel_->text().contains(QStringLiteral("SIMULATOR RESOURCE LIMIT")) ||
        compileStatusLabel_->text().contains(QStringLiteral("OPENED")) || !lastCompilation_.success ||
        circuitSession_ || circuitView_->qubitCount() != 1 || circuitView_->instructionCount() != 2 ||
        !circuitMeasurementBanner_->isHidden() || circuitStateTable_->rowCount() != 0 ||
        circuitStepSlider_->isEnabled() || circuitPlayBtn_->isEnabled() ||
        latestResult_ || !plotResults_.empty() || !currentSweepResults_.empty() ||
        !convergencePlot_->points().empty() || energyGeometry_->text() != QStringLiteral("—")) {
      finishSmokeTest(false, QStringLiteral("simulator initialization failure while opening a bundle left stale output or reported success"));
      return;
    }
    if (!openExperimentFile(smokeQcomp)) {
      finishSmokeTest(false, QStringLiteral("valid saved bundle did not reopen after initialization failure"));
      return;
    }

    // Opening a valid circuit-only bundle clears saved results from the prior H2 run.
    if (!openExperimentFile(initFailurePath) || latestResult_ || !plotResults_.empty() ||
        !currentSweepResults_.empty() || !convergencePlot_->points().empty() ||
        energyGeometry_->text() != QStringLiteral("—") ||
        !progressLabel_->text().contains(QStringLiteral("No chemistry result"))) {
      finishSmokeTest(false, QStringLiteral("circuit-only bundle retained chemistry or convergence data from another run"));
      return;
    }
    QFile::remove(invalidOpenPath);
    QFile::remove(initFailurePath);
    // 1b. Test HeH+ 4-qubit circuit preset in Circuit Workbench
    optimizeCheck_->setChecked(false);
    exampleCombo_->setCurrentIndex(4);
    if (!dslEditor_->toPlainText().contains(QStringLiteral("qubit q0, q1, q2, q3"))) {
      finishSmokeTest(false, QStringLiteral("HeH+ preset did not load 4-qubit declaration"));
      return;
    }
    compileDsl();
    if (!lastCompilation_.success || lastCompilation_.qubitCount != 4) {
      finishSmokeTest(false, QStringLiteral("HeH+ 4-qubit preset did not compile to 4 qubits"));
      return;
    }
    if (circuitSession_->totalSteps() < 50) {
      finishSmokeTest(false, QStringLiteral("HeH+ trial sequence has unexpected instruction count"));
      return;
    }
    const auto layoutMatchesDag = [](const CircuitView &view,
                                     const std::vector<Instruction> &sequence,
                                     int qubits, std::size_t expectedDepth) {
      if (view.instructionCount() != sequence.size() ||
          view.logicalDepth() != expectedDepth) return false;
      const auto dag = compiler::CircuitDag::fromInstructions(
          sequence, static_cast<std::size_t>(qubits));
      if (dag.criticalPathDepth() != expectedDepth) return false;
      for (const compiler::DagNode &node : dag.nodes()) {
        if (node.id >= sequence.size()) return false;
        const auto &layout = view.instructionLayout(node.id);
        if (layout.instructionIndex != node.id ||
            layout.simulationStep != node.id ||
            layout.layerIndex != static_cast<std::size_t>(node.layer)) return false;
      }
      return true;
    };

    compiler::CompilerOptions layoutOptions;
    layoutOptions.optimize = false;
    const auto layoutProbe = compiler::compileProgram(
        "qubit q0, q1, q2, q3\n"
        "CNOT q0, q3\n"
        "X q1\n"
        "H q0\n"
        "X q2\n",
        layoutOptions);
    compiler::CompilerOptions optimizedLayoutOptions;
    optimizedLayoutOptions.optimize = true;
    const auto optimizedLayoutProbe = compiler::compileProgram(
        "qubit q0, q1, q2, q3\n"
        "CNOT q0, q3\n"
        "X q1\n"
        "H q0\n"
        "X q2\n",
        optimizedLayoutOptions);
    CircuitView layoutProbeView;
    layoutProbeView.setCircuit(layoutProbe.rawInstructions, 4);
    if (!layoutProbe.success || !optimizedLayoutProbe.success ||
        !layoutMatchesDag(layoutProbeView, layoutProbe.rawInstructions, 4,
                          layoutProbe.rawDepth) ||
        layoutProbe.rawInstructions.size() != 4 ||
        layoutProbeView.instructionLayout(0).layerIndex != 0 ||
        layoutProbeView.instructionLayout(1).layerIndex != 0 ||
        layoutProbeView.instructionLayout(3).layerIndex != 0 ||
        layoutProbeView.instructionLayout(2).layerIndex != 1 ||
        layoutProbeView.instructionLayout(0).laneIndex ==
            layoutProbeView.instructionLayout(1).laneIndex ||
        layoutProbeView.operationBounds(0).top() > 42 ||
        layoutProbeView.operationBounds(0).bottom() < 42 + 3 * 52) {
      finishSmokeTest(false, QStringLiteral("logical layer layout lost parallel gates, dependencies, or the full CNOT span"));
      return;
    }
    layoutProbeView.setCircuit(optimizedLayoutProbe.instructions, 4);
    if (!layoutMatchesDag(layoutProbeView, optimizedLayoutProbe.instructions, 4,
                          optimizedLayoutProbe.optimizedDepth)) {
      finishSmokeTest(false, QStringLiteral("optimized instructions were not mapped to their own DAG layers and steps"));
      return;
    }

    viewRawRadio_->setChecked(true);
    if (!layoutMatchesDag(*circuitView_, lastCompilation_.rawInstructions, 4,
                          lastCompilation_.rawDepth)) {
      finishSmokeTest(false, QStringLiteral("HeH+ raw display did not map to its sequential instruction indices"));
      return;
    }
    std::size_t widestCnotStep = lastCompilation_.rawInstructions.size();
    int widestCnotSpan = -1;
    for (std::size_t i = 0; i < lastCompilation_.rawInstructions.size(); ++i) {
      const auto &instruction = lastCompilation_.rawInstructions[i];
      if (instruction.gateName == "CNOT" && instruction.qubitIndices.size() >= 2) {
        const int span = std::abs(instruction.qubitIndices[0] -
                                  instruction.qubitIndices[1]);
        if (span > widestCnotSpan) {
          widestCnotSpan = span;
          widestCnotStep = i;
        }
      }
    }
    if (widestCnotStep == lastCompilation_.rawInstructions.size() ||
        widestCnotSpan <= 0) {
      finishSmokeTest(false, QStringLiteral("HeH+ circuit has no inspectable multi-qubit operation"));
      return;
    }
    circuitStepSlider_->setValue(static_cast<int>(widestCnotStep));
    const auto activeOperationVisible = [this](std::size_t step) {
      if (!circuitScrollArea_ || step >= circuitView_->instructionCount() ||
          circuitSession_->currentStep() != step ||
          circuitView_->currentStep() != step ||
          circuitStepSlider_->value() != static_cast<int>(step)) return false;
      const QRect bounds = circuitView_->operationBounds(step);
      QWidget *viewport = circuitScrollArea_->viewport();
      const QRect viewportBounds(QPoint(0, 0), viewport->size());
      const QRect mapped(circuitView_->mapTo(viewport, bounds.topLeft()), bounds.size());
      if (bounds.width() <= viewport->width() &&
          bounds.height() <= viewport->height()) {
        return viewportBounds.contains(mapped);
      }
      return !circuitViewportHint_->isHidden() &&
             viewportBounds.contains(circuitView_->mapTo(viewport, bounds.center()));
    };
    if (circuitSession_->currentStep() != widestCnotStep ||
        !activeOperationVisible(widestCnotStep)) {
      finishSmokeTest(false, QStringLiteral("seeking to a long-circuit CNOT did not expose its complete control/target span"));
      return;
    }
    circuitStepFwdBtn_->click();
    if (!activeOperationVisible(widestCnotStep + 1)) {
      finishSmokeTest(false, QStringLiteral("forward stepping lost active-operation visibility or sequential mapping"));
      return;
    }
    circuitStepBackBtn_->click();
    if (!activeOperationVisible(widestCnotStep)) {
      finishSmokeTest(false, QStringLiteral("backward stepping lost active-operation visibility or sequential mapping"));
      return;
    }

    const auto captureCircuitAt = [this, &activeOperationVisible](int width, int height) {
      showNormal();
      resize(width, height);
      if (mainTabs_) mainTabs_->setCurrentIndex(0);
      QApplication::processEvents();
      const QSize actualSize = size();
      const std::vector<QWidget *> toolbarControls{
          exampleCombo_, compileBtn_, optimizeCheck_, saveExpBtn_, openExpBtn_,
          exportSvgBtn_, exportPngBtn_};
      QWidget *toolbarParent = toolbarControls.front()->parentWidget();
      bool controlsUsable = toolbarParent != nullptr;
      for (std::size_t i = 0; controlsUsable && i < toolbarControls.size(); ++i) {
        const QWidget *control = toolbarControls[i];
        controlsUsable = control && control->isVisible() &&
            control->parentWidget() == toolbarParent &&
            toolbarParent->rect().contains(control->geometry()) &&
            control->width() > 0 && control->height() > 0;
        for (std::size_t j = 0; controlsUsable && j < i; ++j) {
          controlsUsable = !control->geometry().intersects(toolbarControls[j]->geometry());
        }
      }
      if (actualSize != QSize(width, height) || !controlsUsable) {
        std::cout << "QCOMPILE_LAB_WINDOW_GEOMETRY requested="
                  << width << "x" << height << " actual="
                  << actualSize.width() << "x" << actualSize.height()
                  << " controlsUsable=" << controlsUsable << std::endl;
        return false;
      }
      std::cout << "QCOMPILE_LAB_WINDOW_GEOMETRY actual="
                << actualSize.width() << "x" << actualSize.height() << std::endl;
      followCircuitOperation(circuitSession_->currentStep());
      QApplication::processEvents();
      if (!activeOperationVisible(circuitSession_->currentStep())) return false;
      const QString path = QDir::tempPath() +
          QStringLiteral("/qcompile-ms14-circuit-%1-%2.png")
              .arg(width).arg(QDateTime::currentMSecsSinceEpoch());
      if (!grab().save(path)) return false;
      std::cout << "QCOMPILE_MS14_CIRCUIT_SCREENSHOT "
                << path.toStdString() << std::endl;
      return true;
    };
    if (!captureCircuitAt(1024, 680) || !captureCircuitAt(1440, 900)) {
      finishSmokeTest(false, QStringLiteral("main circuit workbench did not fit the requested target size or toolbar controls overlapped"));
      return;
    }

    const int playStart = static_cast<int>(circuitSession_->totalSteps()) - 5;
    circuitStepSlider_->setValue(playStart);
    if (circuitSession_->currentStep() != static_cast<std::size_t>(playStart)) {
      finishSmokeTest(false, QStringLiteral("long-circuit seek did not set the sequential execution step"));
      return;
    }
    circuitPlayBtn_->click();
    QEventLoop playbackWait;
    QTimer::singleShot(650, &playbackWait, &QEventLoop::quit);
    playbackWait.exec();
    circuitPauseBtn_->click();
    if (circuitSession_->currentStep() <= static_cast<std::size_t>(playStart) ||
        !activeOperationVisible(circuitSession_->currentStep())) {
      finishSmokeTest(false, QStringLiteral("playback did not advance and follow the active HeH+ operation"));
      return;
    }
    circuitResetBtn_->click();

    // Step forward 2 instructions: X q0, X q2
    circuitStepFwdBtn_->click();
    circuitStepFwdBtn_->click();
    if (circuitSession_->currentStep() != 2) {
      finishSmokeTest(false, QStringLiteral("HeH+ circuit did not advance 2 steps"));
      return;
    }
    // Index 5 is |0101> in bit order (q0=1, q1=0, q2=1, q3=0)
    if (std::abs(circuitSession_->currentState()[5] - std::complex<double>(1.0, 0.0)) > 1e-12) {
      finishSmokeTest(false, QStringLiteral("HeH+ state after X0, X2 is not determinant |0101>"));
      return;
    }
    // Verify optimization pass reduces gate count on HeH+
    optimizeCheck_->setChecked(true);
    compileDsl();
    if (lastCompilation_.instructions.size() >= lastCompilation_.rawInstructions.size()) {
      finishSmokeTest(false, QStringLiteral("optimization pass did not reduce HeH+ gate count"));
      return;
    }
    viewOptimizedRadio_->setChecked(true);
    if (!layoutMatchesDag(*circuitView_, lastCompilation_.instructions, 4,
                          lastCompilation_.optimizedDepth) ||
        !compileStatusLabel_->text().contains(QStringLiteral("optimized"))) {
      finishSmokeTest(false, QStringLiteral("optimized HeH+ display did not report and render its own logical depth mapping"));
      return;
    }
    viewRawRadio_->setChecked(true);
    if (!layoutMatchesDag(*circuitView_, lastCompilation_.rawInstructions, 4,
                          lastCompilation_.rawDepth)) {
      finishSmokeTest(false, QStringLiteral("switching back to raw HeH+ did not restore its instruction mapping"));
      return;
    }
    viewOptimizedRadio_->setChecked(true);
    circuitResetBtn_->click();
    exampleCombo_->setCurrentIndex(0);

    // 2. Switch to Molecular Chemistry Tab and run existing chemistry smoke verification
    if (mainTabs_) {
      mainTabs_->setCurrentIndex(1);
    }
    if (pesViewCombo_ && energyPlot_) {
      pesViewCombo_->setCurrentIndex(0); // Multi-Molecule Overlay
      if (energyPlot_->viewMode() != 0 || energyPlot_->seriesList().size() < 2) {
        finishSmokeTest(false, QStringLiteral("multi-molecule PES overlay series not properly initialized"));
        return;
      }
      pesViewCombo_->setCurrentIndex(1); // H2 scan
    }

    if (mainTabs_ && comparisonListTable_ && comparisonDeltaTable_) {

    // Comparison identity regression: separate files share the same title and basename.
    QTemporaryDir comparisonRoot;
    const QString firstDirectory = comparisonRoot.filePath(QStringLiteral("first"));
    const QString secondDirectory = comparisonRoot.filePath(QStringLiteral("second"));
    if (!comparisonRoot.isValid() ||
        !QDir().mkpath(firstDirectory) || !QDir().mkpath(secondDirectory)) {
      finishSmokeTest(false, QStringLiteral("could not create comparison identity fixtures"));
      return;
    }
    const QString firstPath = QDir(firstDirectory).filePath(QStringLiteral("run.qcomp"));
    const QString secondPath = QDir(secondDirectory).filePath(QStringLiteral("run.qcomp"));
    try {
      ExperimentBundle firstBundle = loadExperimentBundle(smokeQcomp.toStdString());
      ExperimentBundle secondBundle = firstBundle;
      firstBundle.title = "Shared title";
      secondBundle.title = "Shared title";
      firstBundle.id = "copied-bundle-id";
      secondBundle.id = "copied-bundle-id";
      saveExperimentBundle(firstPath.toStdString(), firstBundle);
      saveExperimentBundle(secondPath.toStdString(), secondBundle);
    } catch (const std::exception &error) {
      finishSmokeTest(false, QString::fromLocal8Bit(error.what()));
      return;
    }

    comparisonBundles_.clear();
    comparisonListTable_->setRowCount(0);
    comparisonDeltaTable_->setRowCount(0);
    comparisonDeltaTable_->horizontalHeaderItem(1)->setText(QStringLiteral("Experiment A"));
    comparisonDeltaTable_->horizontalHeaderItem(2)->setText(QStringLiteral("Experiment B"));
    addExperimentToComparison(firstPath);
    addExperimentToComparison(secondPath);
    const QString firstIdentity = QFileInfo(firstPath).canonicalFilePath();
    const QString secondIdentity = QFileInfo(secondPath).canonicalFilePath();
    if (comparisonListTable_->rowCount() != 2 ||
        comparisonDeltaTable_->rowCount() == 0 ||
        QFileInfo(firstPath).fileName() != QFileInfo(secondPath).fileName() ||
        firstIdentity == secondIdentity ||
        !comparisonListTable_->item(0, 0)->text().contains(QStringLiteral("Shared title")) ||
        !comparisonListTable_->item(1, 0)->text().contains(QStringLiteral("Shared title")) ||
        !comparisonListTable_->item(0, 0)->text().contains(firstIdentity) ||
        !comparisonListTable_->item(1, 0)->text().contains(secondIdentity)) {
      finishSmokeTest(false, QStringLiteral("same-title bundle files were not distinct in the comparison list"));
      return;
    }
    auto *selection = comparisonListTable_->selectionModel();
    selection->clearSelection();
    selection->select(comparisonListTable_->model()->index(0, 0),
                      QItemSelectionModel::Select | QItemSelectionModel::Rows);
    selection->select(comparisonListTable_->model()->index(1, 0),
                      QItemSelectionModel::Select | QItemSelectionModel::Rows);
    const QString resultHeaderA = comparisonDeltaTable_->horizontalHeaderItem(1)->text();
    const QString resultHeaderB = comparisonDeltaTable_->horizontalHeaderItem(2)->text();
    if (selection->selectedRows().size() != 2 ||
        !resultHeaderA.contains(firstIdentity) ||
        !resultHeaderB.contains(secondIdentity) ||
        !comparisonSummaryLabel_->text().contains(firstIdentity) ||
        !comparisonSummaryLabel_->text().contains(secondIdentity)) {
      finishSmokeTest(false, QStringLiteral("comparison selection did not carry distinct file identities into the A/B results"));
      return;
    }
    QFile::remove(smokeQcomp);
    mainTabs_->setCurrentIndex(1);
    }
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
    if (QFileInfo::exists(screenshot)) {
      QFile::remove(screenshot);
    }
    if (width() < 1200 || height() < 800) {
      resize(1440, 920);
    }
    if (!grab().save(screenshot)) success = false;
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

  // Milestone 3: Multi-molecule PES & Live VQE Convergence
  QComboBox *pesViewCombo_ = nullptr;
  ConvergencePlot *convergencePlot_ = nullptr;
  QPushButton *liveH2VqeBtn_ = nullptr;
  QPushButton *liveHeHVqeBtn_ = nullptr;
  QPushButton *cancelVqeBtn_ = nullptr;
  QLabel *vqeStatusLabel_ = nullptr;
  // Interactive Quantum Circuit Workbench widgets & state
  QTabWidget *mainTabs_ = nullptr;
  QPlainTextEdit *dslEditor_ = nullptr;
  QComboBox *exampleCombo_ = nullptr;
  QPushButton *compileBtn_ = nullptr;
  QCheckBox *optimizeCheck_ = nullptr;
  QLabel *compileStatusLabel_ = nullptr;
  QTableWidget *diagnosticsTable_ = nullptr;
  CircuitView *circuitView_ = nullptr;
  QScrollArea *circuitScrollArea_ = nullptr;
  QLabel *circuitViewportHint_ = nullptr;
  QPushButton *circuitResetBtn_ = nullptr;
  QPushButton *circuitStepBackBtn_ = nullptr;
  QPushButton *circuitStepFwdBtn_ = nullptr;
  QPushButton *circuitPlayBtn_ = nullptr;
  QPushButton *circuitPauseBtn_ = nullptr;
  QSlider *circuitStepSlider_ = nullptr;
  QLabel *circuitStepLabel_ = nullptr;
  QRadioButton *viewOptimizedRadio_ = nullptr;
  QRadioButton *viewRawRadio_ = nullptr;
  QLabel *circuitStageLabel_ = nullptr;
  QLabel *circuitMetricsLabel_ = nullptr;
  QLabel *circuitMeasurementBanner_ = nullptr;
  QTableWidget *circuitStateTable_ = nullptr;
  QTimer circuitPlayTimer_;
  QTimer dslDebounceTimer_;
  compiler::CompilationResult lastCompilation_;
  std::unique_ptr<SimulationSession> circuitSession_;
  std::size_t circuitAvailableStepLimit_ = 0;
  bool injectNextCircuitStepResourceFailureForSmoke_ = false;
  bool injectNextSessionInitializationFailureForSmoke_ = false;
  bool isCompiling_ = false;

  // Milestone 4: Persistence, Comparison & Export
  QPushButton *saveExpBtn_ = nullptr;
  QPushButton *openExpBtn_ = nullptr;
  QPushButton *exportSvgBtn_ = nullptr;
  QPushButton *exportPngBtn_ = nullptr;
  QPushButton *exportCsvBtn_ = nullptr;
  QTableWidget *comparisonListTable_ = nullptr;
  QTableWidget *comparisonDeltaTable_ = nullptr;
  QLabel *comparisonSummaryLabel_ = nullptr;
  std::vector<ComparisonEntry> comparisonBundles_;
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
