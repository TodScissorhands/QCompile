#pragma once

#include "compiler/code_generator.h"

#include <QRect>
#include <QWidget>
#include <QString>
#include <cstddef>
#include <vector>

class QPainter;

struct CircuitInstructionLayout {
  std::size_t instructionIndex = 0;
  std::size_t simulationStep = 0;
  std::size_t layerIndex = 0;
  std::size_t laneIndex = 0;
};

class CircuitView : public QWidget {

public:
  explicit CircuitView(QWidget *parent = nullptr);

  void setCircuit(const std::vector<Instruction> &instructions, int qubitCount);
  void setCurrentStep(std::size_t simulationStep);

  int qubitCount() const { return qubitCount_; }
  std::size_t currentStep() const { return currentStep_; }
  std::size_t instructionCount() const { return instructions_.size(); }
  std::size_t logicalDepth() const { return logicalDepth_; }
  const CircuitInstructionLayout &instructionLayout(std::size_t index) const {
    return instructionLayouts_.at(index);
  }
  QRect operationBounds(std::size_t simulationStep) const;

  bool exportPng(const QString &filePath, int scaleFactor = 2);
  void paintCircuitContent(QPainter &painter, int fullWidth, int fullHeight);

  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;
protected:
  void paintEvent(QPaintEvent *event) override;

private:
  void rebuildLayerLayout();
  std::size_t visualColumn(std::size_t instructionIndex) const;
  int wireY(int qubitIndex) const;
  int columnX(std::size_t columnIndex) const;

  std::vector<Instruction> instructions_;
  std::vector<CircuitInstructionLayout> instructionLayouts_;
  std::vector<std::size_t> simulationStepToInstruction_;
  std::vector<std::size_t> layerStartColumns_;
  std::vector<std::size_t> layerColumnCounts_;
  int qubitCount_ = 0;
  std::size_t currentStep_ = 0;
  std::size_t logicalDepth_ = 0;
  std::size_t visualColumnCount_ = 0;
};
