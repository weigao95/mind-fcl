//
// Created by wei on 24-3-30.
//

#pragma once

namespace fcl {

template <typename S>
inline void initializeFixedRotationBoxDisjoint(
    const Transform3<S>& tf_2in1, FixedRotationBoxDisjoint<S>& disjoint) {
  disjoint.rotation_2in1 = tf_2in1.linear().matrix();
  disjoint.translation_2in1 = tf_2in1.translation();
  disjoint.rotation_2in1_abs = disjoint.rotation_2in1.cwiseAbs();
  constexpr S reps = S(1e-6);
  disjoint.rotation_2in1_abs.array() += reps;
}

template <typename S>
inline void FixedRotationBoxDisjoint<S>::initialize(
    const Transform3<S>& tf_2in1) {
  initializeFixedRotationBoxDisjoint(tf_2in1, *this);
}

template <typename S>
inline void FixedRotationBoxDisjoint<S>::initialize(const Transform3<S>& tf_1,
                                                    const Transform3<S>& tf_2) {
  const Transform3<S> tf_2in1 = tf_1.inverse(Eigen::Isometry) * tf_2;
  initialize(tf_2in1);
}

template <typename S>
inline bool FixedRotationBoxDisjoint<S>::isDisjoint(
    const AABB<S>& aabb1, const AABB<S>& aabb2,
    bool check_strict_disjoint) const {
  const Vector3<S> b2_center = aabb2.center();
  const Vector3<S> b2_extent = S(0.5) * (aabb2.max_ - aabb2.min_);
  const Vector3<S> b2_in_b1 =
      rotation_2in1 * b2_center + translation_2in1 - aabb1.center();
  const Vector3<S> b1_extent = S(0.5) * (aabb1.max_ - aabb1.min_);
  return isDisjoint(b2_in_b1, b1_extent, b2_extent, check_strict_disjoint);
}

template <typename S>
bool isOrientedBoxDisjointDefault(const FixedRotationBoxDisjoint<S>& disjoint,
                                  const Vector3<S>& translation_box_2in1,
                                  const Vector3<S>& box1_extent,
                                  const Vector3<S>& box2_extent,
                                  bool check_strict_disjoint) {
  const auto& B = disjoint.rotation_2in1;
  const auto& T = translation_box_2in1;
  const auto& a = box1_extent;
  const auto& b = box2_extent;
  const auto& Bf = disjoint.rotation_2in1_abs;
  S t, s;

  {
    // A1 x A2 = A0
    t = ((T[0] < 0.0) ? -T[0] : T[0]);
    if (t > (a[0] + Bf.row(0).dot(b))) return true;
  }

  {
    // B1 x B2 = B0
    s = B.col(0).dot(T);
    t = ((s < 0.0) ? -s : s);
    if (t > (b[0] + Bf.col(0).dot(a))) return true;
  }

  {
    // A2 x A0 = A1
    t = ((T[1] < 0.0) ? -T[1] : T[1]);
    if (t > (a[1] + Bf.row(1).dot(b))) return true;
  }

  {
    // A0 x A1 = A2
    t = ((T[2] < 0.0) ? -T[2] : T[2]);
    if (t > (a[2] + Bf.row(2).dot(b))) return true;
  }

  {
    // B2 x B0 = B1
    s = B.col(1).dot(T);
    t = ((s < 0.0) ? -s : s);
    if (t > (b[1] + Bf.col(1).dot(a))) return true;
  }

  {
    // B0 x B1 = B2
    s = B.col(2).dot(T);
    t = ((s < 0.0) ? -s : s);
    if (t > (b[2] + Bf.col(2).dot(a))) return true;
  }

  // Not strict
  if (!check_strict_disjoint) return false;

  {
    // A0 x B0
    s = T[2] * B(1, 0) - T[1] * B(2, 0);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[1] * Bf(2, 0) + a[2] * Bf(1, 0) + b[1] * Bf(0, 2) + b[2] * Bf(0, 1)))
      return true;
  }

  {
    // A0 x B1
    s = T[2] * B(1, 1) - T[1] * B(2, 1);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[1] * Bf(2, 1) + a[2] * Bf(1, 1) + b[0] * Bf(0, 2) + b[2] * Bf(0, 0)))
      return true;
  }

  {
    // A0 x B2
    s = T[2] * B(1, 2) - T[1] * B(2, 2);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[1] * Bf(2, 2) + a[2] * Bf(1, 2) + b[0] * Bf(0, 1) + b[1] * Bf(0, 0)))
      return true;
  }

  {
    // A1 x B0
    s = T[0] * B(2, 0) - T[2] * B(0, 0);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(2, 0) + a[2] * Bf(0, 0) + b[1] * Bf(1, 2) + b[2] * Bf(1, 1)))
      return true;
  }

  {
    // A1 x B1
    s = T[0] * B(2, 1) - T[2] * B(0, 1);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(2, 1) + a[2] * Bf(0, 1) + b[0] * Bf(1, 2) + b[2] * Bf(1, 0)))
      return true;
  }

  {
    // A1 x B2
    s = T[0] * B(2, 2) - T[2] * B(0, 2);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(2, 2) + a[2] * Bf(0, 2) + b[0] * Bf(1, 1) + b[1] * Bf(1, 0)))
      return true;
  }

  {
    // A2 x B0
    s = T[1] * B(0, 0) - T[0] * B(1, 0);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(1, 0) + a[1] * Bf(0, 0) + b[1] * Bf(2, 2) + b[2] * Bf(2, 1)))
      return true;
  }

  {
    // A2 x B1
    s = T[1] * B(0, 1) - T[0] * B(1, 1);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(1, 1) + a[1] * Bf(0, 1) + b[0] * Bf(2, 2) + b[2] * Bf(2, 0)))
      return true;
  }

  {
    // A2 x B2
    s = T[1] * B(0, 2) - T[0] * B(1, 2);
    t = ((s < 0.0) ? -s : s);
    if (t >
        (a[0] * Bf(1, 2) + a[1] * Bf(0, 2) + b[0] * Bf(2, 1) + b[1] * Bf(2, 0)))
      return true;
  }

  // Intersect
  return false;
}

template <typename S>
inline bool FixedRotationBoxDisjoint<S>::isDisjoint(
    const Vector3<S>& translation_box_2in1, const Vector3<S>& box1_extent,
    const Vector3<S>& box2_extent, bool check_strict_disjoint) const {
  return isOrientedBoxDisjointDefault(*this, translation_box_2in1, box1_extent,
                                      box2_extent, check_strict_disjoint);
}

// Specialization for float
template <>
inline void FixedRotationBoxDisjoint<float>::initialize(
    const Transform3<float>& tf_2in1) {
  // Init for others
  initializeFixedRotationBoxDisjoint<float>(tf_2in1, *this);

  // Init for simd data
#ifdef FCL_SSE_ENABLED
  const auto& R_r0 = rotation_2in1.row(0);
  rotation_2in1_sse[0] = _mm_setr_ps(R_r0[0], R_r0[1], R_r0[2], 0.f);
  const auto& R_r1 = rotation_2in1.row(1);
  rotation_2in1_sse[1] = _mm_setr_ps(R_r1[0], R_r1[1], R_r1[2], 0.f);
  const auto& R_r2 = rotation_2in1.row(2);
  rotation_2in1_sse[2] = _mm_setr_ps(R_r2[0], R_r2[1], R_r2[2], 0.f);

  constexpr float reps = 1e-6;
  const __m128 epsilonxyz = _mm_setr_ps(reps, reps, reps, 0.f);
  rotation_2in1_abs_sse[0] =
      _mm_add_ps(detail::abs_ps(rotation_2in1_sse[0]), epsilonxyz);
  rotation_2in1_abs_sse[1] =
      _mm_add_ps(detail::abs_ps(rotation_2in1_sse[1]), epsilonxyz);
  rotation_2in1_abs_sse[2] =
      _mm_add_ps(detail::abs_ps(rotation_2in1_sse[2]), epsilonxyz);
#endif
}

#ifdef FCL_SSE_ENABLED

inline bool obbDisjointSSEFloatImpl(const __m128* R, const __m128* AbsR,
                                    const __m128& t, const __m128& r1,
                                    const __m128& r2,
                                    bool check_strict_disjoint) {
  using namespace detail;
  __m128 ra;           // projection of OBB A's halfwidth along three axes
  __m128 rb;           // projection of OBB B's halfwidth along three axes
  __m128 center_dist;  // projection of center distance along three axes

  // Test the three major axes of this OBB.
  if (any_gt_ps(abs_ps(t), _mm_add_ps(r1, mat3x4_mul_vec4(AbsR, r2)))) {
    return true;
  }

  // Test the three major axes of the OBB b.
  center_dist = transp_mat3x4_mul_vec4(R, t);
  if (any_gt_ps(abs_ps(center_dist),
                _mm_add_ps(transp_mat3x4_mul_vec4(AbsR, r1), r2))) {
    return true;
  }

  // Not strict
  if (!check_strict_disjoint) return false;

  // Test the 9 different cross-axes.
  __m128 symmetric_matrix[3] = {
      _mm_setr_ps(0.f, reinterpret_cast<const float*>(&r2)[2],
                  reinterpret_cast<const float*>(&r2)[1], 0.f),
      _mm_setr_ps(reinterpret_cast<const float*>(&r2)[2], 0.f,
                  reinterpret_cast<const float*>(&r2)[0], 0.f),
      _mm_setr_ps(reinterpret_cast<const float*>(&r2)[1],
                  reinterpret_cast<const float*>(&r2)[0], 0.f, 0.f),
  };

  // A.x <cross> B.x
  // A.x <cross> B.y
  // A.x <cross> B.z
  ra = fmadd_ps(vec_splat_ps(r1, 1), AbsR[2],
                _mm_mul_ps(vec_splat_ps(r1, 2), AbsR[1]));
  rb = mat3x4_mul_vec4(symmetric_matrix, AbsR[0]);
  center_dist =
      fmsub_ps(vec_splat_ps(t, 2), R[1], _mm_mul_ps(vec_splat_ps(t, 1), R[2]));
  if (any_gt_ps(abs_ps(center_dist), _mm_add_ps(ra, rb))) {
    return true;
  }

  // A.y <cross> B.x
  // A.y <cross> B.y
  // A.y <cross> B.z
  ra = fmadd_ps(vec_splat_ps(r1, 0), AbsR[2],
                _mm_mul_ps(vec_splat_ps(r1, 2), AbsR[0]));
  rb = mat3x4_mul_vec4(symmetric_matrix, AbsR[1]);
  center_dist =
      fmsub_ps(vec_splat_ps(t, 0), R[2], _mm_mul_ps(vec_splat_ps(t, 2), R[0]));
  if (any_gt_ps(abs_ps(center_dist), _mm_add_ps(ra, rb))) {
    return true;
  }

  // A.z <cross> B.x
  // A.z <cross> B.y
  // A.z <cross> B.z
  ra = fmadd_ps(vec_splat_ps(r1, 0), AbsR[1],
                _mm_mul_ps(vec_splat_ps(r1, 1), AbsR[0]));
  rb = mat3x4_mul_vec4(symmetric_matrix, AbsR[2]);
  center_dist =
      fmsub_ps(vec_splat_ps(t, 1), R[0], _mm_mul_ps(vec_splat_ps(t, 0), R[1]));
  return any_gt_ps(abs_ps(center_dist), _mm_add_ps(ra, rb));
}

#endif

inline bool isOrientedBoxDisjointFloat(
    const FixedRotationBoxDisjoint<float>& disjoint,
    const Vector3<float>& translation_box_2in1,
    const Vector3<float>& box1_extent, const Vector3<float>& box2_extent,
    bool check_strict_disjoint) {
#ifdef FCL_SSE_ENABLED
  const auto& T = translation_box_2in1;
  const auto& a = box1_extent;
  const auto& b = box2_extent;
  __m128 T_sse = _mm_setr_ps(T[0], T[1], T[2], 0.f);
  __m128 a_sse = _mm_setr_ps(a[0], a[1], a[2], 0.f);
  __m128 b_sse = _mm_setr_ps(b[0], b[1], b[2], 0.f);
  return obbDisjointSSEFloatImpl(disjoint.rotation_2in1_sse,
                                 disjoint.rotation_2in1_abs_sse, T_sse, a_sse,
                                 b_sse, check_strict_disjoint);
#else
  return isOrientedBoxDisjointDefault<float>(disjoint, translation_box_2in1,
                                             box1_extent, box2_extent,
                                             check_strict_disjoint);
#endif
}

template <>
inline bool FixedRotationBoxDisjoint<float>::isDisjoint(
    const AABB<float>& aabb1, const AABB<float>& aabb2,
    bool check_strict_disjoint) const {
  // Gather info
  const Vector3f b1_center = aabb1.center();
  const Vector3f b1_extent = 0.5f * (aabb1.max_ - aabb1.min_);
  const Vector3f b2_center = aabb2.center();
  const Vector3f b2_extent = 0.5f * (aabb2.max_ - aabb2.min_);

  // Compute b2_in_b1 =
  //      rotation_2in1 * b2_center + translation_2in1 - aabb1.center();
#ifdef FCL_SSE_ENABLED
  const __m128 b2_center_sse =
      _mm_setr_ps(b2_center[0], b2_center[1], b2_center[2], 0.f);
  const __m128 rotated_b2_center_sse =
      detail::mat3x4_mul_vec4(rotation_2in1_sse, b2_center_sse);
  const float* rotated_b2_center_sse_float =
      reinterpret_cast<const float*>(&rotated_b2_center_sse);
  const float b2_in_b1_x =
      rotated_b2_center_sse_float[0] + translation_2in1[0] - b1_center[0];
  const float b2_in_b1_y =
      rotated_b2_center_sse_float[1] + translation_2in1[1] - b1_center[1];
  const float b2_in_b1_z =
      rotated_b2_center_sse_float[2] + translation_2in1[2] - b1_center[2];
  const Vector3f b2_in_b1(b2_in_b1_x, b2_in_b1_y, b2_in_b1_z);
#else
  const Vector3f b2_in_b1 =
      rotation_2in1 * b2_center + translation_2in1 - b1_center;
#endif
  return isOrientedBoxDisjointFloat(*this, b2_in_b1, b1_extent, b2_extent,
                                    check_strict_disjoint);
}

template <>
inline bool FixedRotationBoxDisjoint<float>::isDisjoint(
    const Vector3<float>& translation_box_2in1,
    const Vector3<float>& box1_extent, const Vector3<float>& box2_extent,
    bool check_strict_disjoint) const {
  return isOrientedBoxDisjointFloat(*this, translation_box_2in1, box1_extent,
                                    box2_extent, check_strict_disjoint);
}

}  // namespace fcl
