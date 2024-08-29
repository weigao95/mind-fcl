//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/common/types.h"
#include "fcl/math/geometry.h"

namespace fcl {

template <typename S>
struct TranslationalDisplacement {
  Vector3<S> unit_axis_in_shape1;
  S scalar_displacement;
};

template <typename S>
struct Interval {
  S lower_bound{-1.0};
  S upper_bound{-1.0};

  // Intersect with another
  void intersect(const Interval& rhs) {
    using std::max;
    using std::min;
    lower_bound = max(lower_bound, rhs.lower_bound);
    upper_bound = min(upper_bound, rhs.upper_bound);
  }

  bool is_empty(S tolerance) const {
    return lower_bound > upper_bound + tolerance;
  }
};

}  // namespace fcl