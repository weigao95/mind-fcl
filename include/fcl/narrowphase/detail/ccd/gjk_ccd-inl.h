#pragma once

namespace fcl {
namespace detail {

template <typename S>
void TranslationalCollisionGJK<S>::setupMinkowskiDiff(
    const GJKGeometryData<S>& swept_volume1, const Transform3<S>& tf1,
    const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
    MinkowskiDiff<S>& shape) {
  // Make mk diff
  shape.shapes[0] = swept_volume1;
  shape.shapes[1] = shape2;
  shape.support_function = detail::computeSupport<S>;
  shape.interior_function = detail::computeInterior<S>;
  shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
  shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;
}

template <typename S>
bool TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionBinary(
    const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
    const Vector3<S>& shape1_displacement_in_tf1,
    const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
    int max_iteration, S tolerance) {
  // Make swept volume
  GJKGeometryData<S> swept_volume1;
  swept_volume1.shape_type = detail::GJKShapeType::ShapeSweptVolume;
  swept_volume1.vec3_data = shape1_displacement_in_tf1;
  swept_volume1.user_ptr = &(shape1);

  // Make mk diff
  MinkowskiDiff<S> shape;
  setupMinkowskiDiff(swept_volume1, tf1, shape2, tf2, shape);

  // Use MPR for binary
  detail::MPR<S> mpr(max_iteration, tolerance);
  auto mpr_result = mpr.Intersect(shape);
  const bool is_overlap_mpr =
      (mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
  return is_overlap_mpr;
}

template <typename S>
bool TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionBinary(
    const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& shape1_displacement,
    const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
    int max_iteration, S tolerance) {
  auto displacement = shape1_displacement.unit_axis_in_shape1 *
                      shape1_displacement.scalar_displacement;
  return CheckSweptVolumeCollisionBinary(shape1, tf1, displacement, shape2, tf2,
                                         max_iteration, tolerance);
}

template <typename S>
std::pair<bool, S>
TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionOneTocSample(
    const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
    const Vector3<S>& shape1_displacement_in_tf1,
    const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
    int max_iteration, S tolerance) {
  // Make swept volume
  GJKGeometryData<S> swept_volume1;
  swept_volume1.shape_type = detail::GJKShapeType::ShapeSweptVolume;
  swept_volume1.vec3_data = shape1_displacement_in_tf1;
  swept_volume1.user_ptr = &(shape1);

  // Make mk diff
  MinkowskiDiff<S> shape;
  setupMinkowskiDiff(swept_volume1, tf1, shape2, tf2, shape);

  // Use MPR for binary
  using MPR = typename detail::MPR<S>;
  using IntersectData = typename MPR::IntersectData;
  using IntersectStatus = typename MPR::IntersectStatus;

  // Run MPR
  MPR mpr(max_iteration, tolerance);
  IntersectData intersect_data;
  IntersectStatus intersect_status = mpr.Intersect(shape, &intersect_data);
  if (intersect_status != IntersectStatus::Intersect) {
    return {false, 0.0};
  }

  // Intersect
  const S toc_sample = computeOneTocSampleForIntersection(
      shape1_displacement_in_tf1, intersect_data);
  return {true, toc_sample};
}

template <typename S>
std::pair<bool, S>
TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionOneTocSample(
    const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& shape1_displacement,
    const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
    int max_iteration, S tolerance) {
  auto displacement = shape1_displacement.unit_axis_in_shape1 *
                      shape1_displacement.scalar_displacement;
  return CheckSweptVolumeCollisionOneTocSample(
      shape1, tf1, displacement, shape2, tf2, max_iteration, tolerance);
}

template <typename S>
S TranslationalCollisionGJK<S>::computeOneTocSampleForIntersection(
    const Vector3<S>& shape1_displacement_in_tf1,
    const typename detail::MPR<S>::IntersectData& intersect_data) {
  const auto& v0 = intersect_data.v0_interior;
  const auto& v1 = intersect_data.v1;
  const auto& v2 = intersect_data.v2;
  const auto& v3 = intersect_data.v3;
  assert(!v0.array().isNaN().any());

  // Degeneration case 1: v0 is the only point
  constexpr S v0_phase = 0.5;
  if (v1.array().isNaN().any()) {
    return v0_phase;
  }

  assert(!v1.array().isNaN().any());
  // Degeneration case 2: one segment of v0 and v1
  if (v2.array().isNaN().any()) {
    // Co-linear and o is in the middle
    const auto& o_to_v0_norm = v0.norm();
    const auto& o_to_v1_norm = v1.norm();
    const auto v0_to_v1_norm = o_to_v0_norm + o_to_v1_norm;

    // Compute weight
    if (v0_to_v1_norm <= 0.0) return v0_phase;
    const auto inv_v0_to_v1_norm = S(1.0) / v0_to_v1_norm;
    const auto v0_weight = o_to_v1_norm * inv_v0_to_v1_norm;
    const auto v1_weight = o_to_v0_norm * inv_v0_to_v1_norm;
    const S v1_phase =
        (intersect_data.v1_dir_in_support.dot(shape1_displacement_in_tf1) > 0)
            ? 1.0
            : 0.0;
    return v0_weight * v0_phase + v1_weight * v1_phase;
  }

  // General case: non NaN
  assert(!v2.array().isNaN().any());
  assert(!v3.array().isNaN().any());

  // Compute the volume
  const auto v0_to_v1 = v1 - v0;
  const auto v0_to_v2 = v2 - v0;
  const auto v0_to_v3 = v3 - v0;
  const S o_v012_volume = std::abs(((v0_to_v1).cross(v0_to_v2)).dot(v0));
  const S o_v013_volume = std::abs(((v0_to_v1).cross(v0_to_v3)).dot(v0));
  const S o_v023_volume = std::abs(((v0_to_v2).cross(v0_to_v3)).dot(v0));
  const S o_v123_volume = std::abs(((v2 - v1).cross(v3 - v1)).dot(v1));
  const S v0123_volume =
      o_v012_volume + o_v013_volume + o_v023_volume + o_v123_volume;

  // Degeneration case?
  if (v0123_volume <= 0.0) {
    // TODO(wei): handle this case
    return 0.0;
  }

  // Volume is positive
  const S inv_v0123_volume = S(1.0) / v0123_volume;

  // Volume should match with total volume
  assert(std::abs(std::abs(((v2 - v1).cross(v3 - v1)).dot(v1 - v0)) -
                  v0123_volume) < 1e-5);

  // Compute the weight
  const S v0_weight = o_v123_volume * inv_v0123_volume;
  const S v1_weight = o_v023_volume * inv_v0123_volume;
  const S v2_weight = o_v013_volume * inv_v0123_volume;
  const S v3_weight = o_v012_volume * inv_v0123_volume;

  // Compute the phase
  const auto& displacement = shape1_displacement_in_tf1;
  const auto& v1_direction = intersect_data.v1_dir_in_support;
  const S v1_phase = (v1_direction.dot(displacement) > 0) ? 1.0 : 0.0;
  const auto& v2_direction = intersect_data.v2_dir_in_support;
  const S v2_phase = (v2_direction.dot(displacement) > 0) ? 1.0 : 0.0;
  const auto& v3_direction = intersect_data.v3_dir_in_support;
  const S v3_phase = (v3_direction.dot(displacement) > 0) ? 1.0 : 0.0;

  // Done
  return v0_weight * v0_phase + v1_weight * v1_phase + v2_weight * v2_phase +
         v3_weight * v3_phase;
}

}  // namespace detail
}  // namespace fcl