//
// Created by Wei Gao on 2024/6/11.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
bool BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
    S box1_half_range_on_axis, S abs_box1_displacement_on_axis,
    bool is_box1_displacement_positive, S box2_half_range_on_axis,
    S box2_offset_axis, Interval<S>& toc) {
  // Using the min/max
  using std::max;
  using std::min;

  // clang-format off
  const S box1_lb = is_box1_displacement_positive ? (-box1_half_range_on_axis) : (-box1_half_range_on_axis - abs_box1_displacement_on_axis);
  const S box1_ub = is_box1_displacement_positive ? (box1_half_range_on_axis + abs_box1_displacement_on_axis) : box1_half_range_on_axis;
  // clang-format on

  // Check bound
  const S box2_lb = box2_offset_axis - box2_half_range_on_axis;
  const S box2_ub = box2_offset_axis + box2_half_range_on_axis;
  if (box1_lb > box2_ub || box2_lb > box1_ub) return true;

  // Compute the interval, which is not empty
  // clang-format off
  if (is_box1_displacement_positive) {
    toc.lower_bound = max(S(0.0), -box1_half_range_on_axis + box2_lb);
    toc.upper_bound = min(abs_box1_displacement_on_axis, box1_half_range_on_axis + box2_ub);
  } else {
    toc.lower_bound = max(S(0.0), -box1_half_range_on_axis - box2_ub);
    toc.upper_bound = min(abs_box1_displacement_on_axis, box1_half_range_on_axis - box2_lb);
  }
  // clang-format on

  // Done
  return false;
}

template <typename S>
bool BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
    const Eigen::Array3<S>& box1_half_range_on_axis,
    const Eigen::Array3<S>& box1_displacement_on_axis,
    const Eigen::Array3<S>& box2_half_range_on_axis,
    const Eigen::Array3<S>& box2_offset_axis, Interval<S>& toc,
    S zero_displacement_tolerance) {
  // clang-format off
  const auto is_box1_displacement_positive = box1_displacement_on_axis > 0.0;
  const auto box1_lb = is_box1_displacement_positive.select(-box1_half_range_on_axis, -box1_half_range_on_axis + box1_displacement_on_axis);
  const auto box1_ub = is_box1_displacement_positive.select(box1_half_range_on_axis + box1_displacement_on_axis, box1_half_range_on_axis);

  // Check bound
  const auto box2_lb = box2_offset_axis - box2_half_range_on_axis;
  const auto box2_ub = box2_offset_axis + box2_half_range_on_axis;
  if ((box1_lb > box2_ub).any() || (box2_lb > box1_ub).any()) {
    return true;
  }

  // Compute toc
  const auto toc_lower_raw = is_box1_displacement_positive.select((-box1_half_range_on_axis + box2_lb), (-box1_half_range_on_axis - box2_ub));
  const auto toc_upper_raw = is_box1_displacement_positive.select((box1_half_range_on_axis + box2_ub), (box1_half_range_on_axis - box2_lb));

  // Scale it
  const auto inv_scale = S(1.0) / box1_displacement_on_axis.abs().max(zero_displacement_tolerance);
  const auto toc_lower_scaled = toc_lower_raw * inv_scale;
  const auto toc_upper_scaled = toc_upper_raw * inv_scale;
  const auto toc_lower = toc_lower_scaled.max(0.0);
  const auto toc_upper = toc_upper_scaled.min(1.0);
  // clang-format on

  // Into toc
  toc.lower_bound = toc_lower.maxCoeff();
  toc.upper_bound = toc_upper.minCoeff();
  return toc.is_empty(interval_empty_tol);
}

template <typename S>
bool BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
    S box1_half_range_on_axis, S box1_displacement_on_axis,
    S box2_half_range_on_axis, S box2_offset_axis, Interval<S>& toc) {
  const bool is_positive_displacement = (box1_displacement_on_axis > 0);
  const auto abs_displacement = is_positive_displacement
                                    ? box1_displacement_on_axis
                                    : -box1_displacement_on_axis;
  return CheckEmptyTocIntervalOnAxis(
      box1_half_range_on_axis, abs_displacement, is_positive_displacement,
      box2_half_range_on_axis, box2_offset_axis, toc);
}

template <typename S>
void BoxPairTranslationalCCD<S>::ScaleIntervalBoxDisjoint(
    Interval<S>& interval_on_axis, S box1_displacement_on_axis_abs,
    S zero_tolerance) {
  if (box1_displacement_on_axis_abs < zero_tolerance) {
    interval_on_axis.lower_bound = 0;
    interval_on_axis.upper_bound = 1;
  } else {
    interval_on_axis.lower_bound /= box1_displacement_on_axis_abs;
    interval_on_axis.upper_bound /= box1_displacement_on_axis_abs;
  }
}

template <typename S>
bool BoxPairTranslationalCCD<S>::isDisjointArray(
    const Vector3<S>& obb1_half_size,
    const TranslationalDisplacement<S>& obb1_displacement,
    const Matrix3<S>& rotation_2in1, const Vector3<S>& translation_2in1,
    const Vector3<S>& obb2_half_size, const Interval<S>& init_toc_interval,
    Interval<S>& interval, S zero_displacement_tolerance) {
  // From obb overlap
  Matrix3<S> rotation_2in1_abs = rotation_2in1.cwiseAbs();
  // const S reps = 1e-6;
  // rotation_2in1_abs.array() += reps;

  const auto scalar_displacement = obb1_displacement.scalar_displacement;
  assert(scalar_displacement >= 0);
  interval = init_toc_interval;

  // Box-1 x/y/z axis
  {
    const auto& box1_half_range_on_axis = obb1_half_size;
    const auto box1_displacement_on_axis =
        obb1_displacement.unit_axis_in_shape1 *
        obb1_displacement.scalar_displacement;
    const auto box2_half_range_on_axis = rotation_2in1_abs * obb2_half_size;
    const auto& box2_offset_axis = translation_2in1;

    // Check disjoint
    Interval<S> toc_interval_box1_xyz;
    const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
        box1_half_range_on_axis.array(), box1_displacement_on_axis.array(),
        box2_half_range_on_axis.array(), box2_offset_axis.array(),
        toc_interval_box1_xyz, zero_displacement_tolerance);
    if (is_empty_interval) {
      return true;
    }

    // Intersect it
    interval.intersect(toc_interval_box1_xyz);
    if (interval.is_empty(interval_empty_tol)) {
      return true;
    }
  }

  // Box-2 x/y/z
  {
    const auto box1_half_range_on_axis =
        rotation_2in1_abs.transpose() * obb1_half_size;
    const auto unit_axis_in_shape2 =
        rotation_2in1.transpose() * obb1_displacement.unit_axis_in_shape1;
    const auto box1_displacement_on_axis =
        unit_axis_in_shape2 * scalar_displacement;
    const auto& box2_half_range_on_axis = obb2_half_size;
    const auto box2_offset_axis = rotation_2in1.transpose() * translation_2in1;

    // Check disjoint
    Interval<S> toc_interval_box2_xyz;
    const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
        box1_half_range_on_axis.array(), box1_displacement_on_axis.array(),
        box2_half_range_on_axis.array(), box2_offset_axis.array(),
        toc_interval_box2_xyz, zero_displacement_tolerance);
    if (is_empty_interval) {
      return true;
    }

    // Intersect it
    interval.intersect(toc_interval_box2_xyz);
    if (interval.is_empty(interval_empty_tol)) {
      return true;
    }
  }

  // Other 9 axis
  for (auto k = 0; k < 3; k++) {
    for (auto i = 0; i < 3; i++) {
      // This axis_i_in_1 is NOT unit, and might be zero actually
      const auto& axis_i_in_1 = Vector3<S>::Unit(k).cross(rotation_2in1.col(i));
      const auto& axis_i_in_1_abs = axis_i_in_1.cwiseAbs();
      const auto box1_half_range_on_axis = axis_i_in_1_abs.dot(obb1_half_size);
      const auto box2_offset_axis = axis_i_in_1.dot(translation_2in1);
      // const auto box1_displacement_on_axis =
      // axis_i_in_1.dot(obb1_displacement.displace_in_shape1);

      // Compute the half range of box2
      const auto& axis_i_in_2 = rotation_2in1.transpose() * axis_i_in_1;
      const auto& axis_i_in_2_abs = axis_i_in_2.cwiseAbs();
      const auto box2_half_range_on_axis = axis_i_in_2_abs.dot(obb2_half_size);

      // Compute displacement on axis
      const auto displacement_axis_projection =
          axis_i_in_1.dot(obb1_displacement.unit_axis_in_shape1);
      const bool is_positive_displacement = displacement_axis_projection > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? displacement_axis_projection * scalar_displacement
              : -displacement_axis_projection * scalar_displacement;

      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
          box1_half_range_on_axis, box1_displacement_on_axis_abs,
          is_positive_displacement, box2_half_range_on_axis, box2_offset_axis,
          toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      ScaleIntervalBoxDisjoint(toc_interval_on_axis,
                               box1_displacement_on_axis_abs,
                               zero_displacement_tolerance);
      interval.intersect(toc_interval_on_axis);
      if (interval.is_empty(interval_empty_tol)) {
        return true;
      }
    }
  }

  // Not disjoint
  assert(!interval.is_empty(interval_empty_tol));
  return false;
}

template <typename S>
bool BoxPairTranslationalCCD<S>::isDisjointNaive(
    const Vector3<S>& obb1_half_size,
    const TranslationalDisplacement<S>& obb1_displacement,
    const Matrix3<S>& rotation_2in1, const Vector3<S>& translation_2in1,
    const Vector3<S>& obb2_half_size, const Interval<S>& init_toc_interval,
    Interval<S>& interval, S zero_displacement_tolerance) {
  // From obb overlap
  Matrix3<S> rotation_2in1_abs = rotation_2in1.cwiseAbs();
  // const S reps = 1e-6;
  // rotation_2in1_abs.array() += reps;

  const auto scalar_displacement = obb1_displacement.scalar_displacement;
  assert(scalar_displacement >= 0);
  interval = init_toc_interval;

  // Box-1 x/y/z axis
  {
    const auto& box1_half_range_on_axis = obb1_half_size;
    // const auto box1_displacement_on_axis =
    // obb1_displacement.unit_axis_in_shape1 *
    // obb1_displacement.scalar_displacement;
    const auto box2_half_range_on_axis = rotation_2in1_abs * obb2_half_size;
    const auto& box2_offset_axis = translation_2in1;
    for (auto i = 0; i < 3; i++) {
      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_positive_displacement =
          obb1_displacement.unit_axis_in_shape1[i] > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? obb1_displacement.unit_axis_in_shape1[i] * scalar_displacement
              : -obb1_displacement.unit_axis_in_shape1[i] * scalar_displacement;
      const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
          box1_half_range_on_axis[i], box1_displacement_on_axis_abs,
          is_positive_displacement, box2_half_range_on_axis[i],
          box2_offset_axis[i], toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      ScaleIntervalBoxDisjoint(toc_interval_on_axis,
                               box1_displacement_on_axis_abs,
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
    const auto box2_offset_axis = rotation_2in1.transpose() * translation_2in1;
    for (auto i = 0; i < 3; i++) {
      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_positive_displacement = unit_axis_in_shape2[i] > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? unit_axis_in_shape2[i] * scalar_displacement
              : -unit_axis_in_shape2[i] * scalar_displacement;
      const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
          box1_half_range_on_axis[i], box1_displacement_on_axis_abs,
          is_positive_displacement, box2_half_range_on_axis[i],
          box2_offset_axis[i], toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      ScaleIntervalBoxDisjoint(toc_interval_on_axis,
                               box1_displacement_on_axis_abs,
                               zero_displacement_tolerance);
      interval.intersect(toc_interval_on_axis);
      if (interval.is_empty(interval_empty_tol)) {
        return true;
      }
    }
  }

  // Other 9 axis
  for (auto k = 0; k < 3; k++) {
    for (auto i = 0; i < 3; i++) {
      // This axis_i_in_1 is NOT unit, and might be zero actually
      const auto& axis_i_in_1 = Vector3<S>::Unit(k).cross(rotation_2in1.col(i));
      const auto& axis_i_in_1_abs = axis_i_in_1.cwiseAbs();
      const auto box1_half_range_on_axis = axis_i_in_1_abs.dot(obb1_half_size);
      const auto box2_offset_axis = axis_i_in_1.dot(translation_2in1);
      // const auto box1_displacement_on_axis =
      // axis_i_in_1.dot(obb1_displacement.displace_in_shape1);

      // Compute the half range of box2
      const auto& axis_i_in_2 = rotation_2in1.transpose() * axis_i_in_1;
      const auto& axis_i_in_2_abs = axis_i_in_2.cwiseAbs();
      const auto box2_half_range_on_axis = axis_i_in_2_abs.dot(obb2_half_size);

      // Compute displacement on axis
      const auto displacement_axis_projection =
          axis_i_in_1.dot(obb1_displacement.unit_axis_in_shape1);
      const bool is_positive_displacement = displacement_axis_projection > 0;
      const S box1_displacement_on_axis_abs =
          is_positive_displacement
              ? displacement_axis_projection * scalar_displacement
              : -displacement_axis_projection * scalar_displacement;

      // Compute interval
      Interval<S> toc_interval_on_axis;
      const bool is_empty_interval = CheckEmptyTocIntervalOnAxis(
          box1_half_range_on_axis, box1_displacement_on_axis_abs,
          is_positive_displacement, box2_half_range_on_axis, box2_offset_axis,
          toc_interval_on_axis);
      if (is_empty_interval) {
        return true;
      }

      // Scale the interval
      ScaleIntervalBoxDisjoint(toc_interval_on_axis,
                               box1_displacement_on_axis_abs,
                               zero_displacement_tolerance);
      interval.intersect(toc_interval_on_axis);
      if (interval.is_empty(interval_empty_tol)) {
        return true;
      }
    }
  }

  // Not disjoint
  assert(!interval.is_empty(interval_empty_tol));
  return false;
}

template <typename S>
bool BoxPairTranslationalCCD<S>::IsDisjoint(
    const OBB<S>& obb1, const TranslationalDisplacement<S>& obb1_displacement,
    const OBB<S>& obb2, Interval<S>& interval, S zero_displacement_tolerance) {
  Interval<S> init_interval;
  init_interval.lower_bound = 0.0;
  init_interval.upper_bound = 1.0;
  return IsDisjoint(obb1, obb1_displacement, obb2, init_interval, interval,
                    zero_displacement_tolerance);
}

template <typename S>
bool BoxPairTranslationalCCD<S>::IsDisjoint(
    const OBB<S>& obb1, const TranslationalDisplacement<S>& obb1_displacement,
    const OBB<S>& obb2, const Interval<S>& init_toc_interval_bound,
    Interval<S>& interval, S zero_displacement_tolerance) {
  const auto& box1_rotation = obb1.axis;
  const auto& box2_rotation = obb2.axis;
  const Vector3<S> t_in_world = obb2.To - obb1.To;
  const auto t_in_box1 = box1_rotation.transpose() * t_in_world;
  return isDisjointArray(obb1.extent, obb1_displacement,
                         box1_rotation.transpose() * box2_rotation, t_in_box1,
                         obb2.extent, init_toc_interval_bound, interval,
                         zero_displacement_tolerance);
}

}  // namespace detail
}  // namespace fcl
