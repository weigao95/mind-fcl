#pragma once
#include "fcl/narrowphase/detail/primitive_shape_algorithm/simplex_common.h"

namespace fcl {

namespace detail {

template <typename S>
inline bool isBoxTriangleSeparatedByTriangleAxis(
    const Vector3<S>& box_half_size, const Vector3<S>& triangle_edge,
    const Vector3<S>& p_on_edge, const Vector3<S>& p_not_on_edge) {
  const auto& e0 = triangle_edge;
  const auto& v0 = p_on_edge;
  const auto& v1 = p_not_on_edge;
  const S e0_x = e0[0];
  const S e0_y = e0[1];
  const S e0_z = e0[2];
  using std::abs;
  const S abs_e0_x = abs(e0_x);
  const S abs_e0_y = abs(e0_y);
  const S abs_e0_z = abs(e0_z);

  // The axis0 is cross(box_x, e0), [0, e0_z, -e0_y]
  {
    const S box_radius_axis =
        box_half_size[1] * abs_e0_z + box_half_size[2] * abs_e0_y;
    const auto v0_dot_axis = v0[1] * e0_z - v0[2] * e0_y;
    const auto v1_dot_axis = v1[1] * e0_z - v1[2] * e0_y;
    if (v0_dot_axis < v1_dot_axis) {
      if (v0_dot_axis > box_radius_axis) return true;
      if (v1_dot_axis < -box_radius_axis) return true;
    } else {
      if (v1_dot_axis > box_radius_axis) return true;
      if (v0_dot_axis < -box_radius_axis) return true;
    }
  }

  // The axis1 is cross(box_y, e0), [e0_z, 0, -e0_x]
  {
    const S box_radius_axis =
        box_half_size[0] * abs_e0_z + box_half_size[2] * abs_e0_x;
    const auto v0_dot_axis = v0[0] * e0_z - v0[2] * e0_x;
    const auto v1_dot_axis = v1[0] * e0_z - v1[2] * e0_x;
    if (v0_dot_axis < v1_dot_axis) {
      if (v0_dot_axis > box_radius_axis) return true;
      if (v1_dot_axis < -box_radius_axis) return true;
    } else {
      if (v1_dot_axis > box_radius_axis) return true;
      if (v0_dot_axis < -box_radius_axis) return true;
    }
  }

  // The axis1 is cross(box_z, e0), [e0_y, -e0_x, 0]
  {
    const S box_radius_axis =
        box_half_size[0] * abs_e0_y + box_half_size[1] * abs_e0_x;
    const auto v0_dot_axis = v0[0] * e0_y - v0[1] * e0_x;
    const auto v1_dot_axis = v1[0] * e0_y - v1[1] * e0_x;
    if (v0_dot_axis < v1_dot_axis) {
      if (v0_dot_axis > box_radius_axis) return true;
      if (v1_dot_axis < -box_radius_axis) return true;
    } else {
      if (v1_dot_axis > box_radius_axis) return true;
      if (v0_dot_axis < -box_radius_axis) return true;
    }
  }

  // Cannot determine separated
  return false;
}

template <typename S>
bool detectBoxTriangleOverlap(const Vector3<S>& box_half_size,
                              const Vector3<S>& v0, const Vector3<S>& v1,
                              const Vector3<S>& v2) {
  // Test on box axis x/y/z
  {
    S min_v, max_v;
    find_min_max(v0[0], v1[0], v2[0], min_v, max_v);
    if (max_v < -box_half_size[0] || min_v > box_half_size[0]) return false;
    find_min_max(v0[1], v1[1], v2[1], min_v, max_v);
    if (max_v < -box_half_size[1] || min_v > box_half_size[1]) return false;
    find_min_max(v0[2], v1[2], v2[2], min_v, max_v);
    if (max_v < -box_half_size[2] || min_v > box_half_size[2]) return false;
  }

  // Edge tests
  const Vector3<S> e0 = v1 - v0;
  const bool separated_by_e0 =
      isBoxTriangleSeparatedByTriangleAxis<S>(box_half_size, e0, v0, v2);
  if (separated_by_e0) return false;

  const Vector3<S> e1 = v2 - v1;
  const bool separated_by_e1 =
      isBoxTriangleSeparatedByTriangleAxis<S>(box_half_size, e1, v1, v0);
  if (separated_by_e1) return false;

  const Vector3<S> e2 = v0 - v2;
  const bool separated_by_e2 =
      isBoxTriangleSeparatedByTriangleAxis<S>(box_half_size, e2, v0, v1);
  if (separated_by_e2) return false;

  // face test
  const Vector3<S> triangle_normal = e0.cross(e1);
  Vector3<S> p_min;
  Vector3<S> p_max;
  for (int i = 0; i < 3; ++i) {
    const auto point = v0[i];
    if (triangle_normal[i] > S(0.0)) {
      p_min[i] = -box_half_size[i] - point;
      p_max[i] = box_half_size[i] - point;
    } else {
      p_min[i] = box_half_size[i] - point;
      p_max[i] = -box_half_size[i] - point;
    }
  }
  return triangle_normal.dot(p_min) <= 0 && triangle_normal.dot(p_max) >= 0;
}

template <typename S>
bool boxTriangleIntersectImpl(const Box<S>& s, const Vector3<S>& P1,
                              const Vector3<S>& P2, const Vector3<S>& P3,
                              const Transform3<S>& tf_triangle2box) {
  const auto v1 = tf_triangle2box * P1;
  const auto v2 = tf_triangle2box * P2;
  const auto v3 = tf_triangle2box * P3;
  const auto box_half_size = S(0.5) * s.side;
  return detectBoxTriangleOverlap<S>(box_half_size, v1, v2, v3);
}

template <typename S>
bool boxTriangleIntersect(const Box<S>& s, const Transform3<S>& tf1,
                          const Vector3<S>& P1, const Vector3<S>& P2,
                          const Vector3<S>& P3, const Transform3<S>& tf2) {
  return boxTriangleIntersectImpl(s, P1, P2, P3,
                                  tf1.inverse(Eigen::Isometry) * tf2);
}

template <typename S>
bool boxTriangleIntersect(const Box<S>& s, const Transform3<S>& tf1,
                          const Vector3<S>& P1, const Vector3<S>& P2,
                          const Vector3<S>& P3) {
  return boxTriangleIntersectImpl(s, P1, P2, P3, tf1.inverse(Eigen::Isometry));
}

}  // namespace detail
}  // namespace fcl
