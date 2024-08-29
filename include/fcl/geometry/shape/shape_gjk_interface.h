//
// Created by wei on 23-3-18.
//

#pragma once

#include "fcl/common/types.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/capsule.h"
#include "fcl/geometry/shape/cone.h"
#include "fcl/geometry/shape/convex.h"
#include "fcl/geometry/shape/cylinder.h"
#include "fcl/geometry/shape/ellipsoid.h"
#include "fcl/geometry/shape/halfspace.h"
#include "fcl/geometry/shape/plane.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/geometry/shape/triangle_p.h"
#include "fcl/geometry/shape/tetrahedron.h"

namespace fcl {
namespace detail {

// Convert fcl::Shape to GJKGeometry
template <typename Shape>
GJKGeometryData<typename Shape::S> shapeToGJK(const Shape& shape);

// New shape type
template <typename T>
GJKGeometryData<T> triangleToGJK(const TriangleP<T>* triangle);
template <typename T>
Vector3<T> triangleSupport(const GJKGeometryData<T>& geometry_data,
                           const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> tetrahedronToGJK(const Tetrahedron<T>* tetrahedron);
template <typename T>
Vector3<T> tetrahedronSupport(const GJKGeometryData<T>& geometry_data,
                              const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> convexToGJK(const Convex<T>* convex);
template <typename T>
Vector3<T> convexSupport(const GJKGeometryData<T>& geometry_data,
                         const Vector3<T>& direction);

// Construct existing types from fcl::ShapeBase
template <typename T>
GJKGeometryData<T> boxToGJK(const Box<T>& box);

template <typename T>
GJKGeometryData<T> sphereToGJK(const Sphere<T>& sphere);

template <typename T>
GJKGeometryData<T> ellipsoidToGJK(const Ellipsoid<T>& ellipsoid);

template <typename T>
GJKGeometryData<T> capsuleToGJK(const Capsule<T>& capsule);

template <typename T>
GJKGeometryData<T> coneToGJK(const Cone<T>& cone);

template <typename T>
GJKGeometryData<T> cylinderToGJK(const Cylinder<T>& cylinder);

}  // namespace detail
}  // namespace fcl

#include "fcl/geometry/shape/shape_gjk_interface-inl.h"
