#include "compiler/coupling_map.h"

#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>

namespace compiler {

CouplingMap::CouplingMap(std::size_t qubitCount)
    : qubitCount_(qubitCount), edgeCount_(0), adjacencyList_(qubitCount) {}

bool CouplingMap::addEdge(std::size_t u, std::size_t v) {
  if (u >= qubitCount_ || v >= qubitCount_ || u == v) {
    return false;
  }

  if (hasEdge(u, v)) {
    return true; // Idempotent duplicate edge addition
  }

  adjacencyList_[u].push_back(v);
  adjacencyList_[v].push_back(u);

  // Preserve deterministic ascending neighbor order for BFS tie-breaking
  std::sort(adjacencyList_[u].begin(), adjacencyList_[u].end());
  std::sort(adjacencyList_[v].begin(), adjacencyList_[v].end());

  ++edgeCount_;
  return true;
}

bool CouplingMap::hasEdge(std::size_t u, std::size_t v) const noexcept {
  if (u >= qubitCount_ || v >= qubitCount_ || u == v) {
    return false;
  }
  const auto &uNeighbors = adjacencyList_[u];
  return std::binary_search(uNeighbors.begin(), uNeighbors.end(), v);
}

const std::vector<std::size_t> &CouplingMap::neighbors(std::size_t u) const {
  if (u >= qubitCount_) {
    throw std::out_of_range("CouplingMap: qubit index out of bounds");
  }
  return adjacencyList_[u];
}

bool CouplingMap::isConnected() const {
  if (qubitCount_ <= 1) {
    return true;
  }

  std::vector<bool> visited(qubitCount_, false);
  std::queue<std::size_t> q;

  visited[0] = true;
  q.push(0);
  std::size_t visitedCount = 0;

  while (!q.empty()) {
    std::size_t curr = q.front();
    q.pop();
    ++visitedCount;

    for (std::size_t neighbor : adjacencyList_[curr]) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        q.push(neighbor);
      }
    }
  }

  return visitedCount == qubitCount_;
}

bool CouplingMap::findShortestPath(
    std::size_t source,
    std::size_t target,
    std::vector<std::size_t> &path
) const {
  path.clear();

  if (source >= qubitCount_ || target >= qubitCount_) {
    return false;
  }

  if (source == target) {
    path.push_back(source);
    return true;
  }

  std::vector<bool> visited(qubitCount_, false);
  std::vector<std::size_t> parent(qubitCount_, 0);
  std::queue<std::size_t> q;

  visited[source] = true;
  q.push(source);

  bool found = false;
  while (!q.empty()) {
    std::size_t curr = q.front();
    q.pop();

    if (curr == target) {
      found = true;
      break;
    }

    // Neighbors are visited in ascending numerical order (deterministic tie-breaking)
    for (std::size_t neighbor : adjacencyList_[curr]) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        parent[neighbor] = curr;
        q.push(neighbor);
      }
    }
  }

  if (!found) {
    return false;
  }

  // Backtrack to reconstruct shortest path from target to source
  std::size_t curr = target;
  while (curr != source) {
    path.push_back(curr);
    curr = parent[curr];
  }
  path.push_back(source);
  std::reverse(path.begin(), path.end());

  return true;
}

CouplingMap CouplingMap::createLine(std::size_t qubitCount) {
  if (qubitCount == 0) {
    throw std::invalid_argument("CouplingMap: physical qubit count must be positive");
  }
  CouplingMap map(qubitCount);
  for (std::size_t i = 0; i + 1 < qubitCount; ++i) {
    map.addEdge(i, i + 1);
  }
  return map;
}

CouplingMap CouplingMap::createRing(std::size_t qubitCount) {
  if (qubitCount == 0) {
    throw std::invalid_argument("CouplingMap: physical qubit count must be positive");
  }
  CouplingMap map = createLine(qubitCount);
  if (qubitCount > 2) {
    map.addEdge(qubitCount - 1, 0);
  }
  return map;
}

CouplingMap CouplingMap::createStar(std::size_t qubitCount, std::size_t center) {
  if (qubitCount == 0) {
    throw std::invalid_argument("CouplingMap: physical qubit count must be positive");
  }
  if (center >= qubitCount) {
    throw std::invalid_argument("CouplingMap: star center out of bounds");
  }
  CouplingMap map(qubitCount);
  for (std::size_t i = 0; i < qubitCount; ++i) {
    if (i != center) {
      map.addEdge(center, i);
    }
  }
  return map;
}

CouplingMap CouplingMap::createGrid(std::size_t rows, std::size_t cols) {
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("CouplingMap: grid dimensions must be positive");
  }
  if (rows > std::numeric_limits<std::size_t>::max() / cols) {
    throw std::overflow_error("CouplingMap: grid size overflows size_t");
  }

  std::size_t qubitCount = rows * cols;
  CouplingMap map(qubitCount);

  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      std::size_t u = r * cols + c;
      if (c + 1 < cols) {
        map.addEdge(u, r * cols + (c + 1));
      }
      if (r + 1 < rows) {
        map.addEdge(u, (r + 1) * cols + c);
      }
    }
  }

  return map;
}

CouplingMap CouplingMap::createFullyConnected(std::size_t qubitCount) {
  if (qubitCount == 0) {
    throw std::invalid_argument("CouplingMap: physical qubit count must be positive");
  }
  CouplingMap map(qubitCount);
  for (std::size_t i = 0; i < qubitCount; ++i) {
    for (std::size_t j = i + 1; j < qubitCount; ++j) {
      map.addEdge(i, j);
    }
  }
  return map;
}

} // namespace compiler
