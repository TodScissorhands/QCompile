#include "gui/circuit_view.h"

#include "compiler/circuit_dag.h"

#include <QFont>
#include <QFontMetrics>
#include <QPaintEvent>
#include <QPainter>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace {

constexpr int kWireSpacing = 52;
constexpr int kLeftMargin = 64;
constexpr int kTopMargin = 42;
constexpr int kColumnSpacing = 68;

struct WireSpan {
  int first = 0;
  int last = 0;
};

WireSpan wireSpan(const Instruction &instruction, int qubitCount) {
  int first = std::max(0, qubitCount);
  int last = -1;
  for (const int qubit : instruction.qubitIndices) {
    if (qubit < 0 || qubit >= qubitCount) continue;
    first = std::min(first, qubit);
    last = std::max(last, qubit);
  }
  if (last < 0) {
    return {0, std::max(0, qubitCount - 1)};
  }
  return {first, last};
}

bool spansOverlap(const WireSpan &first, const WireSpan &second) {
  return first.first <= second.last && second.first <= first.last;
}

QString formatGateLabel(const Instruction &instr) {
  if (instr.gateName == "RZ") {
    std::ostringstream ss;
    ss << "RZ(" << std::fixed << std::setprecision(2) << instr.parameterValue << ")";
    return QString::fromStdString(ss.str());
  }
  return QString::fromStdString(instr.gateName);
}

} // namespace

CircuitView::CircuitView(QWidget *parent) : QWidget(parent) {
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
  setMinimumSize(minimumSizeHint());
}

void CircuitView::setCircuit(const std::vector<Instruction> &instructions, int qubitCount) {
  instructions_ = instructions;
  qubitCount_ = std::max(0, qubitCount);
  currentStep_ = 0;
  rebuildLayerLayout();
  updateGeometry();
  update();
}

void CircuitView::setCurrentStep(std::size_t simulationStep) {
  currentStep_ = std::min(simulationStep, simulationStepToInstruction_.size());
  update();
}

void CircuitView::rebuildLayerLayout() {
  instructionLayouts_.assign(instructions_.size(), {});
  simulationStepToInstruction_.assign(instructions_.size(), instructions_.size());
  layerStartColumns_.clear();
  layerColumnCounts_.clear();
  logicalDepth_ = 0;
  visualColumnCount_ = 0;
  if (instructions_.empty()) return;

  const compiler::CircuitDag dag = compiler::CircuitDag::fromInstructions(
      instructions_, static_cast<std::size_t>(qubitCount_));
  logicalDepth_ = dag.criticalPathDepth();
  std::vector<std::vector<std::vector<WireSpan>>> occupiedSpans(logicalDepth_);

  // fromInstructions() appends one DAG node per input instruction in order;
  // DagNode::id is therefore the index in this exact raw or optimized sequence.
  // Use node identity rather than layerSchedule()'s copied Instructions, which
  // cannot distinguish duplicate gate values at different simulation steps.
  for (const compiler::DagNode &node : dag.nodes()) {
    if (node.id >= instructions_.size() || node.layer < 0 ||
        static_cast<std::size_t>(node.layer) >= logicalDepth_) {
      continue;
    }
    const std::size_t layer = static_cast<std::size_t>(node.layer);
    const WireSpan span = wireSpan(instructions_[node.id], qubitCount_);
    std::size_t lane = 0;
    for (; lane < occupiedSpans[layer].size(); ++lane) {
      const bool overlaps = std::any_of(
          occupiedSpans[layer][lane].begin(), occupiedSpans[layer][lane].end(),
          [&span](const WireSpan &other) { return spansOverlap(span, other); });
      if (!overlaps) break;
    }
    if (lane == occupiedSpans[layer].size()) {
      occupiedSpans[layer].emplace_back();
    }
    occupiedSpans[layer][lane].push_back(span);
    instructionLayouts_[node.id] = {node.id, node.id, layer, lane};
    simulationStepToInstruction_[node.id] = node.id;
  }

  layerStartColumns_.resize(logicalDepth_);
  layerColumnCounts_.resize(logicalDepth_);
  for (std::size_t layer = 0; layer < logicalDepth_; ++layer) {
    layerStartColumns_[layer] = visualColumnCount_;
    layerColumnCounts_[layer] = std::max<std::size_t>(1, occupiedSpans[layer].size());
    visualColumnCount_ += layerColumnCounts_[layer];
  }
}

std::size_t CircuitView::visualColumn(std::size_t instructionIndex) const {
  const CircuitInstructionLayout &layout = instructionLayouts_.at(instructionIndex);
  return layerStartColumns_.at(layout.layerIndex) + layout.laneIndex;
}

int CircuitView::wireY(int qubitIndex) const {
  return kTopMargin + qubitIndex * kWireSpacing;
}

int CircuitView::columnX(std::size_t columnIndex) const {
  return kLeftMargin + static_cast<int>(columnIndex) * kColumnSpacing + 32;
}

QRect CircuitView::operationBounds(std::size_t simulationStep) const {
  const std::size_t instructionIndex =
      simulationStepToInstruction_.at(simulationStep);
  const Instruction &instruction = instructions_.at(instructionIndex);
  const int x = columnX(visualColumn(instructionIndex));
  if (instruction.type == InstructionType::Measure) {
    QRect bounds;
    for (const int qubit : instruction.qubitIndices) {
      if (qubit < 0 || qubit >= qubitCount_) continue;
      const QRect gateBounds(x - 20, wireY(qubit) - 14, 40, 28);
      bounds = bounds.isNull() ? gateBounds : bounds.united(gateBounds);
    }
    if (!bounds.isNull()) return bounds;
  } else if (instruction.gateName == "CNOT" &&
             instruction.qubitIndices.size() >= 2) {
    const int control = instruction.qubitIndices[0];
    const int target = instruction.qubitIndices[1];
    if (control < 0 || control >= qubitCount_ ||
        target < 0 || target >= qubitCount_) return {};
    const int controlY = wireY(control);
    const int targetY = wireY(target);
    return QRect(x - 12, std::min(controlY, targetY) - 12, 24,
                 std::abs(controlY - targetY) + 24);
  } else if (!instruction.qubitIndices.empty()) {
    const int qubit = instruction.qubitIndices[0];
    if (qubit < 0 || qubit >= qubitCount_) return {};
    const QString label = formatGateLabel(instruction);
    const QFontMetrics metrics(QFont(QStringLiteral("Monospace"), 9, QFont::Bold));
    const int boxWidth = std::max(40, metrics.horizontalAdvance(label) + 12);
    return QRect(x - boxWidth / 2, wireY(qubit) - 14, boxWidth, 28);
  }
  return QRect(x - 18, kTopMargin - 14, 36, 28);
}

QSize CircuitView::sizeHint() const {
  const int width = std::max(
      520, kLeftMargin + static_cast<int>(visualColumnCount_ + 1) * kColumnSpacing + 48);
  const int height = std::max(
      180, kTopMargin + std::max(1, qubitCount_) * kWireSpacing + 28);
  return QSize(width, height);
}

QSize CircuitView::minimumSizeHint() const {
  return QSize(460, 160);
}

bool CircuitView::exportPng(const QString &filePath, int scaleFactor) {
  if (qubitCount_ <= 0 || instructions_.empty()) return false;
  const QSize size = sizeHint();
  const int scale = std::max(1, scaleFactor);
  const int width = size.width() * scale;
  const int height = size.height() * scale;

  QImage image(width, height, QImage::Format_ARGB32);
  image.fill(QColor(11, 15, 23));
  {
    QPainter painter(&image);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.scale(scale, scale);
    paintCircuitContent(painter, size.width(), size.height());
  }
  return image.save(filePath);
}

void CircuitView::paintEvent(QPaintEvent *) {
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);
  paintCircuitContent(painter, std::max(width(), sizeHint().width()),
                      std::max(height(), sizeHint().height()));
}

void CircuitView::paintCircuitContent(QPainter &painter, int fullWidth, int fullHeight) {
  painter.fillRect(0, 0, fullWidth, fullHeight, QColor(11, 15, 23));

  if (qubitCount_ <= 0) {
    painter.setPen(QColor(100, 116, 139));
    painter.setFont(QFont(QStringLiteral("Monospace"), 10));
    painter.drawText(QRect(0, 0, fullWidth, fullHeight), Qt::AlignCenter,
                     QStringLiteral("No quantum register declared (compile a circuit to visualize tracks)"));
    return;
  }

  const int diagramBottom = wireY(qubitCount_ - 1) + 22;
  for (std::size_t layer = 0; layer < logicalDepth_; ++layer) {
    const std::size_t firstColumn = layerStartColumns_[layer];
    const std::size_t lastColumn = firstColumn + layerColumnCounts_[layer] - 1;
    const int firstX = columnX(firstColumn);
    const int lastX = columnX(lastColumn);
    const int left = firstX - kColumnSpacing / 2;
    const int right = lastX + kColumnSpacing / 2;
    const QRect layerRect(left, kTopMargin - 20, right - left, diagramBottom - (kTopMargin - 20));

    painter.fillRect(layerRect, layer % 2 == 0
        ? QColor(15, 23, 42, 72) : QColor(30, 41, 59, 40));
    painter.setPen(QPen(QColor(42, 58, 82), 1, Qt::DotLine));
    painter.drawLine(left, layerRect.top(), left, layerRect.bottom());
    painter.setPen(QColor(100, 116, 139));
    painter.setFont(QFont(QStringLiteral("Monospace"), 7, QFont::Bold));
    painter.drawText(QRect(left, 2, right - left, 13), Qt::AlignCenter,
                     QStringLiteral("LAYER %1").arg(layer + 1));
  }
  if (logicalDepth_ > 0) {
    const int right = columnX(visualColumnCount_ - 1) + kColumnSpacing / 2;
    painter.setPen(QPen(QColor(42, 58, 82), 1, Qt::DotLine));
    painter.drawLine(right, kTopMargin - 20, right, diagramBottom);
  }

  // Qubit wires and labels.
  painter.setFont(QFont(QStringLiteral("Monospace"), 9, QFont::Bold));
  for (int qubit = 0; qubit < qubitCount_; ++qubit) {
    const int y = wireY(qubit);
    painter.setPen(QPen(QColor(30, 41, 59), 2, Qt::SolidLine));
    painter.drawLine(kLeftMargin - 12, y, fullWidth - 20, y);

    const QRect labelRect(12, y - 10, 36, 20);
    painter.setPen(QPen(QColor(51, 65, 85), 1));
    painter.setBrush(QColor(15, 23, 42));
    painter.drawRoundedRect(labelRect, 2, 2);
    painter.setPen(QColor(148, 163, 184));
    painter.drawText(labelRect, Qt::AlignCenter, QStringLiteral("q%1").arg(qubit));
  }

  // The cursor is a sequential instruction step, not a logical layer cursor.
  if (currentStep_ < simulationStepToInstruction_.size()) {
    const std::size_t instructionIndex =
        simulationStepToInstruction_.at(currentStep_);
    const int cursorX = columnX(visualColumn(instructionIndex));
    const QRectF highlightRect(cursorX - 26, kTopMargin - 2, 52,
                               diagramBottom - (kTopMargin - 2));
    painter.setPen(QPen(QColor(56, 189, 248, 180), 1.5, Qt::DashLine));
    painter.setBrush(QColor(56, 189, 248, 24));
    painter.drawRoundedRect(highlightRect, 4, 4);

    const QRectF pillRect(cursorX - 22, kTopMargin - 20, 44, 14);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(56, 189, 248));
    painter.drawRoundedRect(pillRect, 2, 2);
    painter.setPen(QColor(11, 15, 23));
    painter.setFont(QFont(QStringLiteral("Monospace"), 7, QFont::Black));
    painter.drawText(pillRect, Qt::AlignCenter, QStringLiteral("NEXT"));
  }

  for (std::size_t index = 0; index < instructions_.size(); ++index) {
    const Instruction &instruction = instructions_[index];
    const CircuitInstructionLayout &layout = instructionLayouts_[index];
    const int x = columnX(visualColumn(index));
    const bool isActive = layout.simulationStep == currentStep_;
    const bool isExecuted = layout.simulationStep < currentStep_;

    QColor borderColor = isExecuted ? QColor(71, 85, 105) : QColor(30, 41, 59);
    QColor fillColor = isExecuted ? QColor(30, 41, 59) : QColor(15, 23, 42);
    QColor textColor = isExecuted ? QColor(241, 245, 249) : QColor(148, 163, 184);
    if (isActive) {
      borderColor = QColor(56, 189, 248);
      fillColor = QColor(15, 23, 42);
      textColor = QColor(56, 189, 248);
    }

    if (instruction.type == InstructionType::Measure) {
      for (const int qubit : instruction.qubitIndices) {
        if (qubit < 0 || qubit >= qubitCount_) continue;
        const int y = wireY(qubit);
        const QRect gateRect(x - 20, y - 14, 40, 28);
        painter.setPen(QPen(isActive ? QColor(129, 140, 248) : QColor(99, 102, 241),
                            isActive ? 2 : 1));
        painter.setBrush(QColor(30, 27, 75));
        painter.drawRoundedRect(gateRect, 4, 4);
        painter.setPen(QColor(199, 210, 254));
        painter.drawText(gateRect, Qt::AlignCenter, QStringLiteral("◷ M"));
      }
    } else if (instruction.gateName == "CNOT" &&
               instruction.qubitIndices.size() >= 2) {
      const int control = instruction.qubitIndices[0];
      const int target = instruction.qubitIndices[1];
      if (control >= 0 && control < qubitCount_ &&
          target >= 0 && target < qubitCount_) {
        const int controlY = wireY(control);
        const int targetY = wireY(target);
        painter.setPen(QPen(borderColor, isActive ? 2.5 : 1.5, Qt::SolidLine));
        painter.drawLine(x, controlY, x, targetY);
        painter.setPen(Qt::NoPen);
        painter.setBrush(borderColor);
        painter.drawEllipse(QPointF(x, controlY), 5.0, 5.0);
        constexpr double targetRadius = 10.0;
        painter.setPen(QPen(borderColor, isActive ? 2 : 1.5));
        painter.setBrush(fillColor);
        painter.drawEllipse(QPointF(x, targetY), targetRadius, targetRadius);
        painter.drawLine(x - targetRadius, targetY, x + targetRadius, targetY);
        painter.drawLine(x, targetY - targetRadius, x, targetY + targetRadius);
      }
    } else if (!instruction.qubitIndices.empty()) {
      const int qubit = instruction.qubitIndices[0];
      if (qubit >= 0 && qubit < qubitCount_) {
        const int y = wireY(qubit);
        const QString label = formatGateLabel(instruction);
        const QFontMetrics metrics(painter.font());
        const int boxWidth = std::max(40, metrics.horizontalAdvance(label) + 12);
        const QRect gateRect(x - boxWidth / 2, y - 14, boxWidth, 28);
        painter.setPen(QPen(borderColor, isActive ? 2 : 1.2));
        painter.setBrush(fillColor);
        painter.drawRoundedRect(gateRect, 4, 4);
        painter.setPen(textColor);
        painter.drawText(gateRect, Qt::AlignCenter, label);
      }
    }
  }
}
