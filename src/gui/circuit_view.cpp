#include "gui/circuit_view.h"

#include <QPainter>
#include <QPaintEvent>
#include <QFont>
#include <QFontMetrics>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace {

constexpr int kWireSpacing = 52;
constexpr int kLeftMargin = 64;
constexpr int kTopMargin = 42;
constexpr int kColumnSpacing = 68;

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
  qubitCount_ = qubitCount;
  currentStep_ = 0;
  updateGeometry();
  update();
}

void CircuitView::setCurrentStep(std::size_t step) {
  currentStep_ = step;
  update();
}

int CircuitView::wireY(int qubitIndex) const {
  return kTopMargin + qubitIndex * kWireSpacing;
}

int CircuitView::columnX(std::size_t columnIndex) const {
  return kLeftMargin + static_cast<int>(columnIndex) * kColumnSpacing + 32;
}

QSize CircuitView::sizeHint() const {
  const int w = std::max(520, kLeftMargin + static_cast<int>(instructions_.size() + 1) * kColumnSpacing + 48);
  const int h = std::max(180, kTopMargin + std::max(1, qubitCount_) * kWireSpacing + 28);
  return QSize(w, h);
}

QSize CircuitView::minimumSizeHint() const {
  return QSize(460, 160);
}

void CircuitView::paintEvent(QPaintEvent *) {
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);

  // Background
  painter.fillRect(rect(), QColor(11, 15, 23));

  if (qubitCount_ <= 0) {
    painter.setPen(QColor(100, 116, 139));
    painter.setFont(QFont(QStringLiteral("Monospace"), 10));
    painter.drawText(rect(), Qt::AlignCenter,
                     QStringLiteral("No quantum register declared (compile a circuit to visualize tracks)"));
    return;
  }

  // 1. Draw horizontal qubit wires and labels
  painter.setFont(QFont(QStringLiteral("Monospace"), 9, QFont::Bold));
  const int fullWidth = std::max(width(), sizeHint().width());

  for (int q = 0; q < qubitCount_; ++q) {
    const int y = wireY(q);

    // Qubit wire track
    painter.setPen(QPen(QColor(30, 41, 59), 2, Qt::SolidLine));
    painter.drawLine(kLeftMargin - 12, y, fullWidth - 20, y);

    // Qubit label badge
    const QRect labelRect(12, y - 10, 36, 20);
    painter.setPen(QPen(QColor(51, 65, 85), 1));
    painter.setBrush(QColor(15, 23, 42));
    painter.drawRoundedRect(labelRect, 2, 2);

    painter.setPen(QColor(148, 163, 184));
    painter.drawText(labelRect, Qt::AlignCenter, QStringLiteral("q%1").arg(q));
  }

  // 2. Active step cursor column highlight
  if (!instructions_.empty() && currentStep_ < instructions_.size()) {
    const int cursorX = columnX(currentStep_);
    const int topY = kTopMargin - 22;
    const int bottomY = wireY(qubitCount_ - 1) + 22;
    const QRectF highlightRect(cursorX - 26, topY, 52, bottomY - topY);

    painter.setPen(QPen(QColor(56, 189, 248, 180), 1.5, Qt::DashLine));
    painter.setBrush(QColor(56, 189, 248, 24));
    painter.drawRoundedRect(highlightRect, 4, 4);

    // Cursor header pill
    const QRectF pillRect(cursorX - 22, topY - 14, 44, 14);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(56, 189, 248));
    painter.drawRoundedRect(pillRect, 2, 2);
    painter.setPen(QColor(11, 15, 23));
    painter.setFont(QFont(QStringLiteral("Monospace"), 7, QFont::Black));
    painter.drawText(pillRect, Qt::AlignCenter, QStringLiteral("NEXT"));
  }

  // 3. Render gate instructions
  painter.setFont(QFont(QStringLiteral("Monospace"), 8, QFont::Bold));

  for (std::size_t i = 0; i < instructions_.size(); ++i) {
    const auto &instr = instructions_[i];
    const int x = columnX(i);
    const bool isActive = (i == currentStep_);
    const bool isExecuted = (i < currentStep_);

    QColor borderColor = isExecuted ? QColor(71, 85, 105) : QColor(30, 41, 59);
    QColor fillColor = isExecuted ? QColor(30, 41, 59) : QColor(15, 23, 42);
    QColor textColor = isExecuted ? QColor(241, 245, 249) : QColor(148, 163, 184);

    if (isActive) {
      borderColor = QColor(56, 189, 248);
      fillColor = QColor(15, 23, 42);
      textColor = QColor(56, 189, 248);
    }

    if (instr.type == InstructionType::Measure) {
      // Measurement instruction
      for (int qubit : instr.qubitIndices) {
        if (qubit < 0 || qubit >= qubitCount_) continue;
        const int y = wireY(qubit);
        const QRect gateRect(x - 20, y - 14, 40, 28);

        painter.setPen(QPen(isActive ? QColor(129, 140, 248) : QColor(99, 102, 241), isActive ? 2 : 1));
        painter.setBrush(QColor(30, 27, 75));
        painter.drawRoundedRect(gateRect, 4, 4);

        painter.setPen(QColor(199, 210, 254));
        painter.drawText(gateRect, Qt::AlignCenter, QStringLiteral("◷ M"));
      }
    } else if (instr.gateName == "CNOT" && instr.qubitIndices.size() >= 2) {
      // CNOT: control dot and target circle with cross
      const int cQ = instr.qubitIndices[0];
      const int tQ = instr.qubitIndices[1];
      if (cQ >= 0 && cQ < qubitCount_ && tQ >= 0 && tQ < qubitCount_) {
        const int cY = wireY(cQ);
        const int tY = wireY(tQ);

        // Vertical connecting bus
        painter.setPen(QPen(borderColor, isActive ? 2.5 : 1.5, Qt::SolidLine));
        painter.drawLine(x, cY, x, tY);

        // Control dot (solid disc)
        painter.setPen(Qt::NoPen);
        painter.setBrush(borderColor);
        painter.drawEllipse(QPointF(x, cY), 5.0, 5.0);

        // Target symbol: circle with plus crosshair
        const double targetRadius = 10.0;
        painter.setPen(QPen(borderColor, isActive ? 2 : 1.5));
        painter.setBrush(fillColor);
        painter.drawEllipse(QPointF(x, tY), targetRadius, targetRadius);

        painter.drawLine(x - targetRadius, tY, x + targetRadius, tY);
        painter.drawLine(x, tY - targetRadius, x, tY + targetRadius);
      }
    } else {
      // Single-qubit unitary gate (H, X, Z, RZ)
      if (!instr.qubitIndices.empty()) {
        const int qubit = instr.qubitIndices[0];
        if (qubit >= 0 && qubit < qubitCount_) {
          const int y = wireY(qubit);
          const QString label = formatGateLabel(instr);
          const QFontMetrics fm(painter.font());
          const int boxW = std::max(40, fm.horizontalAdvance(label) + 12);
          const QRect gateRect(x - boxW / 2, y - 14, boxW, 28);

          painter.setPen(QPen(borderColor, isActive ? 2 : 1.2));
          painter.setBrush(fillColor);
          painter.drawRoundedRect(gateRect, 4, 4);

          painter.setPen(textColor);
          painter.drawText(gateRect, Qt::AlignCenter, label);
        }
      }
    }
  }
}
