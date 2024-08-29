#pragma once

#include "fcl/geometry/shape/box.h"

namespace fcl {

namespace detail {
template <typename S>
bool boxTriangleIntersect(const Box<S>& s, const Transform3<S>& tf1,
                          const Vector3<S>& P1, const Vector3<S>& P2,
                          const Vector3<S>& P3, const Transform3<S>& tf2);

template <typename S>
bool boxTriangleIntersect(const Box<S>& s, const Transform3<S>& tf1,
                          const Vector3<S>& P1, const Vector3<S>& P2,
                          const Vector3<S>& P3);

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_triangle-inl.h"
