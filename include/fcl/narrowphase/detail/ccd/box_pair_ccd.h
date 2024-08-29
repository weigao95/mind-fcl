//
// Created by Wei Gao on 2024/6/11.
//

#pragma once

#include "fcl/math/bv/OBB.h"
#include "fcl/narrowphase/detail/ccd/ccd_typedef.h"

namespace fcl {
namespace detail {

template <typename S>
struct BoxPairTranslationalCCD {
  static bool CheckEmptyTocIntervalOnAxis(S box1_half_range_on_axis,
                                          S box1_displacement_on_axis,
                                          S box2_half_range_on_axis,
                                          S box2_offset_axis, Interval<S>& toc);
  static bool CheckEmptyTocIntervalOnAxis(S box1_half_range_on_axis,
                                          S abs_box1_displacement_on_axis,
                                          bool is_box1_displacement_positive,
                                          S box2_half_range_on_axis,
                                          S box2_offset_axis, Interval<S>& toc);
  static bool CheckEmptyTocIntervalOnAxis(
      const Eigen::Array3<S>& box1_half_range_on_axis,
      const Eigen::Array3<S>& box1_displacement_on_axis,
      const Eigen::Array3<S>& box2_half_range_on_axis,
      const Eigen::Array3<S>& box2_offset_axis, Interval<S>& toc,
      S zero_tolerance);


  static bool IsDisjoint(const OBB<S>& obb1,
                         const TranslationalDisplacement<S>& obb1_displacement,
                         const OBB<S>& obb2, Interval<S>& toc_interval,
                         S zero_displacement_tolerance);
  static bool IsDisjoint(const OBB<S>& obb1,
                         const TranslationalDisplacement<S>& obb1_displacement,
                         const OBB<S>& obb2,
                         const Interval<S>& init_toc_interval_bound,
                         Interval<S>& toc_interval,
                         S zero_displacement_tolerance);
  static void ScaleIntervalBoxDisjoint(Interval<S>& interval_on_axis,
                                       S box1_displacement_on_axis_abs,
                                       S zero_tolerance);

  // Private utility
 private:
  // Constants
  static constexpr S interval_empty_tol = 0.0;

  // Naive impl of obb disjoint
  static bool isDisjointNaive(
      const Vector3<S>& obb1_half_size,
      const TranslationalDisplacement<S>& obb1_displacement,
      const Matrix3<S>& rotation_2in1, const Vector3<S>& translation_2in1,
      const Vector3<S>& obb2_half_size, const Interval<S>& init_toc_interval,
      Interval<S>& toc_interval, S zero_displacement_tolerance);
  static bool isDisjointArray(
      const Vector3<S>& obb1_half_size,
      const TranslationalDisplacement<S>& obb1_displacement,
      const Matrix3<S>& rotation_2in1, const Vector3<S>& translation_2in1,
      const Vector3<S>& obb2_half_size, const Interval<S>& init_toc_interval,
      Interval<S>& toc_interval, S zero_displacement_tolerance);
};

template <typename S>
struct FixedOrientationBoxPairTranslationalCCD {
  // Basic meta
  TranslationalDisplacement<S> obb1_displacement;
  Matrix3<S> rotation_2in1;
  Vector3<S> frame_translation_2in1;
  Vector3<S> unit_axis_in_box2;
  Matrix3<S> rotation_2in1_abs;

  // Axis in 1 and 2
  static constexpr int kNumCrossAxis = 9;
  std::array<Vector3<S>, kNumCrossAxis> cross_axis_in_1;
  // std::array<Vector3<S>, kNumCrossAxis> cross_axis_in_2;
  std::array<Vector3<S>, kNumCrossAxis> abs_cross_axis_in_1;
  std::array<Vector3<S>, kNumCrossAxis> abs_cross_axis_in_2;
  std::array<S, kNumCrossAxis> unit_displacement_dot_cross_axis;

  inline void Initialize(
      const Transform3<S>& tf_2in1,
      const TranslationalDisplacement<S>& obb1_displacement_in);
  inline void Initialize(
      const Transform3<S>& tf_1, const Transform3<S>& tf_2,
      const TranslationalDisplacement<S>& obb1_displacement_in);
  inline bool IsDisjoint(const AABB<S>& aabb1, const AABB<S>& aabb2,
                         Interval<S>& interval,
                         S zero_displacement_tolerance) const;
  inline bool IsDisjoint(const AABB<S>& aabb1, const AABB<S>& aabb2,
                         const Interval<S>& init_interval,
                         Interval<S>& interval,
                         S zero_displacement_tolerance) const;
  inline bool IsDisjoint(const Vector3<S>& obb1_half_size,
                         const Vector3<S>& translation_box_2in1,
                         const Vector3<S>& obb2_half_size,
                         const Interval<S>& init_interval,
                         Interval<S>& interval,
                         S zero_displacement_tolerance) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/box_pair_ccd-inl.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd_fixed_orientation-inl.h"