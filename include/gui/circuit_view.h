#pragma once

#include "compiler/code_generator.h"

#include <QWidget>
#include <QString>
#include <cstddef>
#include <vector>

class QPainter;

class CircuitView : public QWidget {

public:
  explicit CircuitView(QWidget *parent = nullptr);

  void setCircuit(const std::vector<Instruction> &instructions, int qubitCount);
  void setCurrentStep(std::size_t step);

  int qubitCount() const { return qubitCount_; }
  std::size_t currentStep() const { return currentStep_; }
  std::size_t instructionCount() const { return instructions_.size(); }
  bool exportPng(const QString &filePath, int scaleFactor = 2);
  void paintCircuitContent(QPainter &painter, int fullWidth, int fullHeight);

  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;
protected:
  void paintEvent(QPaintEvent *event) override;

private:
  std::vector<Instruction> instructions_;
  int qubitCount_ = 0;
  std::size_t currentStep_ = 0;

  int wireY(int qubitIndex) const;
  int columnX(std::size_t columnIndex) const;
};
