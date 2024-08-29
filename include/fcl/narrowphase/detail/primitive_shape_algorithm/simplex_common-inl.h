#pragma once

namespace fcl {
namespace detail {

template <typename S>
inline void find_min_max(S x0, S x1, S x2, S& min, S& max) {
  min = max = x0;

  if (x1 < min) {
    min = x1;
  }
  if (x1 > max) {
    max = x1;
  }
  if (x2 < min) {
    min = x2;
  }
  if (x2 > max) {
    max = x2;
  }
}

template <typename S>
inline void find_min_max(S x0, S x1, S x2, S x3, S& min, S& max) {
  min = max = x0;

  if (x1 < min) {
    min = x1;
  }
  if (x1 > max) {
    max = x1;
  }
  if (x2 < min) {
    min = x2;
  }
  if (x2 > max) {
    max = x2;
  }
  if (x3 < min) {
    min = x3;
  }
  if (x3 > max) {
    max = x3;
  }
}

}  // namespace detail
}  // namespace fcl