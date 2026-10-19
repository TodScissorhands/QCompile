#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace compiler {

// Byte offsets and columns use UTF-8 bytes. Lines and columns are 1-based;
// offsets are 0-based. A SourceSpan covers [start, end), so an insertion point
// has identical start and end positions. A zero line/column is unlocated.
struct SourcePosition {
  std::size_t offset = 0;
  int line = 0;
  int column = 0;
};

struct SourceSpan {
  SourcePosition start;
  SourcePosition end;

  bool isValid() const noexcept {
    return start.line > 0 && start.column > 0 && end.line > 0 &&
           end.column > 0 && end.offset >= start.offset &&
           (end.line > start.line ||
            (end.line == start.line && end.column >= start.column));
  }
};

inline SourceSpan pointSpan(SourcePosition position) noexcept {
  return {position, position};
}

inline SourceSpan joinedSpan(const SourceSpan &first,
                             const SourceSpan &last) noexcept {
  if (!first.isValid() || !last.isValid()) return {};
  return {first.start, last.end};
}

class SourceLocatedError : public std::runtime_error {
public:
  SourceLocatedError(std::string message, SourceSpan span)
      : std::runtime_error(std::move(message)), span_(span) {}

  const SourceSpan &span() const noexcept { return span_; }

private:
  SourceSpan span_;
};

} // namespace compiler

using compiler::SourceLocatedError;
using compiler::SourcePosition;
using compiler::SourceSpan;
using compiler::joinedSpan;
using compiler::pointSpan;
