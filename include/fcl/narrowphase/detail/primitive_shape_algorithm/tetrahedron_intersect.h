#pragma once

#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/tetrahedron.h"
#include "fcl/geometry/shape/triangle_p.h"

namespace fcl {

namespace detail {
template <typename S>
bool tetrahedronTrahedronIntersect(const Tetrahedron<S>& s1,
                                   const Transform3<S>& tf1,
                                   const Tetrahedron<S>& s2,
                                   const Transform3<S>& tf2);

template <typename S>
bool boxTerahedronIntersect(const Box<S>& s1, const Transform3<S>& tf1,
                           const Tetrahedron<S>& s2, const Transform3<S>& tf2);

template <typename S>
bool triangleTerahedronIntersect(const TriangleP<S>& s1,
                                const Transform3<S>& tf1,
                                const Tetrahedron<S>& s2,
                                const Transform3<S>& tf2);
}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/primitive_shape_algorithm/tetrahedron_intersect-inl.h"
