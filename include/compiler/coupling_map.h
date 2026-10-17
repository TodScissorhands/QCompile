#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace compiler {

class CouplingMap {
public:
  CouplingMap() = default;
  explicit CouplingMap(std::size_t qubitCount);

  std::size_t qubitCount() const noexcept { return qubitCount_; }
  std::size_t edgeCount() const noexcept { return edgeCount_; }

  bool addEdge(std::size_t u, std::size_t v);
  bool hasEdge(std::size_t u, std::size_t v) const noexcept;
  const std::vector<std::size_t> &neighbors(std::size_t u) const;
  bool isConnected() const;

  // Shortest path via BFS with deterministic ascending neighbor tie-breaking
  bool findShortestPath(
      std::size_t source,
      std::size_t target,
      std::vector<std::size_t> &path
  ) const;

  // Validated Topology Factories
  static CouplingMap createLine(std::size_t qubitCount);
  static CouplingMap createRing(std::size_t qubitCount);
  static CouplingMap createStar(std::size_t qubitCount, std::size_t center = 0);
  static CouplingMap createGrid(std::size_t rows, std::size_t cols);
  static CouplingMap createFullyConnected(std::size_t qubitCount);

private:
  std::size_t qubitCount_ = 0;
  std::size_t edgeCount_ = 0;
  std::vector<std::vector<std::size_t>> adjacencyList_;
};

} // namespace compiler
