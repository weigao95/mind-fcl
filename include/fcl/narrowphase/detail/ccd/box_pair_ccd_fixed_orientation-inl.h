//
// Created by Wei Gao on 2024/6/14.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
void FixedOrientationBoxPairTranslationalCCD<S>::Initialize(
    const Transform3<S>& tf_1, const Transform3<S>& tf_2,
    const TranslationalDisplacement<S>& obb1_displacement_in) {
  const Transform3<S> tf_2in1 = tf_1.inverse(Eigen::Isometry) * tf_2;
  Initialize(tf_2in1, obb1_displacement_in);
}

template <typename S>
void FixedOrientationBoxPairTranslationalCCD<S>::Initialize(
    const Transform3<S>& tf_2in1,
    const TranslationalDisplacement<S>& obb1_displacement_in) {
  obb1_displacement = obb1_displacement_in;
  rotation_2in1 = tf_2in1.linear().matrix();
  frame_translation_2in1 = tf_2in1.translation();
  const auto& unit_displacement_dir_in_1 =
      obb1_displacement.unit_axis_in_shape1;
  unit_axis_in_box2 = rotation_2in1 * unit_displacement_dir_in_1;

  // Abs rotation for range computation
  rotation_2in1_abs = rotation_2in1.cwiseAbs();
  // constexpr S reps = S(1e-6);
  // rotation_2in1_abs.array() += reps;

  // Compute the axis
  int offset = 0;
  for (auto k = 0; k < 3; k++) {
    for (auto i = 0; i < 3; i++) {
      const auto new_axis_in_1 =
          Vector3<S>::Unit(k).cross(rotation_2in1.col(i));
      cross_axis_in_1[offset] = new_axis_in_1;
      unit_displacement_dot_cross_axis[offset] =
          new_axis_in_1.dot(obb1_displacement.unit_axis_in_shape1);
      abs_cross_axis_in_1[offset] = new_axis_in_1.cwiseAbs();

      // Info frame 2
      const auto& new_axis_in_2 = rotation_2in1.transpose() * new_axis_in_1;
      // cross_axis_in_2[offset] = new_axis_in_2;
      abs_cross_axis_in_2[offset] = new_axis_in_2.cwiseAbs();

      // Update offset
      offset += 1;
    }
  }
}

template <typename S>
bool FixedOrientationBoxPairTranslationalCCD<S>::IsDisjoint(
    const AABB<S>& aabb1, const AABB<S>& aabb2, Interval<S>& interval,
    S zero_displacement_tolerance) const {
  // Into interface
  Interval<S> init_interval;
  init_interval.lower_bound = 0.0;
  init_interval.upper_bound = 1.0;
  return IsDisjoint(aabb1, aabb2, init_interval, interval,
                    zero_displacement_tolerance);
}

template <typename S>
bool FixedOrientationBoxPairTranslationalCCD<S>::IsDisjoint(
    const AABB<S>& aabb1, const AABB<S>& aabb2,
    const Interval<S>& init_interval, Interval<S>& interval,
    S zero_displacement_tolerance) const {
  // Gather info
  const Vector3<S> b1_center = aabb1.center();
  const Vector3<S> b1_half_size = S(0.5) * (aabb1.max_ - aabb1.min_);
  const Vector3<S> b2_center = aabb2.center();
  const Vector3<S> b2_half_size = S(0.5) * (aabb2.max_ - aabb2.min_);
  const Vector3<S> b2_in_b1 =
      rotation_2in1 * b2_center + frame_translation_2in1 - b1_center;

  // Into interface
  return IsDisjoint(b1_half_size, b2_in_b1, b2_half_size, init_interval,
                    interval, zero_displacement_tolerance);
}

template <typename S>
bool FixedOrientationBoxPairTranslationalCCD<S>::IsDisjoint(
    const Vector3<S>& obb1_half_size, const Vector3<S>& translation_box_2in1,
    const Vector3<S>& obb2_half_size, const Interval<S>& init_interval,
    Interval<S>& interval, S zero_displacement_tolerance) const {
  const auto scalar_displacement = obb1_displacement.scalar_displacement;
  assert(scalar_displacement >= 0);
  interval = init_interval;
  constexpr S interval_empty_tol = 0.0;

  // Box-1 x/y/z axis
  {
    const auto& box1_half_range_on_axis = obb1_half_size;
    // const auto box1_displacement_on_axis =
    // obb1_displacement.unit_axis_in_shape1 *
    // obb1_displacement.scalar_displacement;
    const auto box2_half_range_on_axis = rotation_2in1_abs * obb2_half_size;
    const auto& box2_offset_axis = translation_box_2in1;
    for (auto i = 0; i < 3; i++) {
      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_positive_displacement =
          obb1_displacement.unit_axis_in_shape1[i] > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? obb1_displacement.unit_axis_in_shape1[i] * scalar_displacement
              : -obb1_displacement.unit_axis_in_shape1[i] * scalar_displacement;
      const bool is_empty_interval =
          BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
              box1_half_range_on_axis[i], box1_displacement_on_axis_abs,
              is_positive_displacement, box2_half_range_on_axis[i],
              box2_offset_axis[i], toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      BoxPairTranslationalCCD<S>::ScaleIntervalBoxDisjoint(
          toc_interval_on_axis, box1_displacement_on_axis_abs,
          zero_displacement_tolerance);
      interval.intersect(toc_interval_on_axis);
      if (interval.is_empty(interval_empty_tol)) {
        return true;
      }
    }
  }

  // Box-2 x/y/z
  {
    const auto box1_half_range_on_axis =
        rotation_2in1_abs.transpose() * obb1_half_size;
    const auto unit_axis_in_shape2 =
        rotation_2in1.transpose() * obb1_displacement.unit_axis_in_shape1;
    // const auto box1_displacement_on_axis = unit_axis_in_shape2 *
    // scalar_displacement;
    const auto& box2_half_range_on_axis = obb2_half_size;
    const auto box2_offset_axis =
        rotation_2in1.transpose() * translation_box_2in1;
    for (auto i = 0; i < 3; i++) {
      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_positive_displacement = unit_axis_in_shape2[i] > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? unit_axis_in_shape2[i] * scalar_displacement
              : -unit_axis_in_shape2[i] * scalar_displacement;
      const bool is_empty_interval =
          BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
              box1_half_range_on_axis[i], box1_displacement_on_axis_abs,
              is_positive_displacement, box2_half_range_on_axis[i],
              box2_offset_axis[i], toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      BoxPairTranslationalCCD<S>::ScaleIntervalBoxDisjoint(
          toc_interval_on_axis, box1_displacement_on_axis_abs,
          zero_displacement_tolerance);
      interval.intersect(toc_interval_on_axis);
      if (interval.is_empty(interval_empty_tol)) {
        return true;
      }
    }
  }

  // Other 9 axis
  for (auto i = 0; i < kNumCrossAxis; i++) {
    // Gather data
    const auto& axis_i_in_1 = cross_axis_in_1[i];
    const auto& axis_i_in_1_abs = abs_cross_axis_in_1[i];
    const auto& axis_i_in_2_abs = abs_cross_axis_in_2[i];
    const auto displacement_axis_projection =
        unit_displacement_dot_cross_axis[i];

    // Compute range
    const auto box1_half_range_on_axis = axis_i_in_1_abs.dot(obb1_half_size);
    const auto box2_offset_axis = axis_i_in_1.dot(translation_box_2in1);
    const auto box2_half_range_on_axis = axis_i_in_2_abs.dot(obb2_half_size);

    // Compute displacement on axis
    const bool is_positive_displacement = displacement_axis_projection > 0;
    const S box1_displacement_on_axis_abs =
        is_positive_displacement
            ? displacement_axis_projection * scalar_displacement
            : -displacement_axis_projection * scalar_displacement;

    // Compute interval
    Interval<S> toc_interval_on_axis;
    const bool is_empty_interval =
        BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
            box1_half_range_on_axis, box1_displacement_on_axis_abs,
            is_positive_displacement, box2_half_range_on_axis, box2_offset_axis,
            toc_interval_on_axis);
    if (is_empty_interval) {
      return true;
    }

    // Scale the interval
    BoxPairTranslationalCCD<S>::ScaleIntervalBoxDisjoint(
        toc_interval_on_axis, box1_displacement_on_axis_abs,
        zero_displacement_tolerance);
    interval.intersect(toc_interval_on_axis);
    if (interval.is_empty(interval_empty_tol)) {
      return true;
    }
  }

  // Not disjoint
  assert(!interval.is_empty(interval_empty_tol));
  return false;
}

}  // namespace detail
}  // namespace fcl