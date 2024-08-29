//
// Created by wei on 24-3-30.
//

#pragma once

#include "fcl/math/bv/AABB.h"
#include "fcl/math/math_simd_details.h"

namespace fcl {

template <typename S>
struct FixedRotationBoxDisjoint {
  // Construct once, and test multiple times
  inline void initialize(const Transform3<S>& tf_2in1);
  inline void initialize(const Transform3<S>& tf_1, const Transform3<S>& tf_2);
  inline bool isDisjoint(const Vector3<S>& translation_box_2in1,
                  const Vector3<S>& box1_extent, const Vector3<S>& box2_extent,
                  bool check_strict_disjoint) const;
  inline bool isDisjoint(const AABB<S>& aabb1, const AABB<S>& aabb2,
                  bool check_strict_disjoint) const;

  // Raw data
  Matrix3<S> rotation_2in1;
  Vector3<S> translation_2in1;
  Matrix3<S> rotation_2in1_abs;

#ifdef FCL_SSE_ENABLED
  __m128 rotation_2in1_sse[3]; // For float
  __m128 rotation_2in1_abs_sse[3];
#endif
};

}  // namespace fcl

#include "fcl/math/fixed_rotation_obb_disjoint-inl.h"
