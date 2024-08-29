//
// Created by Wei Gao on 2023/3/18.
//

#pragma once

#include "fcl/geometry/shape/shape_gjk_interface.h"

namespace fcl {
namespace detail {

/// The customized interface for fcl::ShapeBase
template <typename T>
GJKGeometryData<T> constructGJKGeometry(const ShapeBase<T>* shape);
template <typename T>
Vector3<T> computeSupport(const GJKGeometryData<T>& gjk_geometry,
                          const Vector3<T>& direction);
template <typename T>
Vector3<T> computeInterior(const GJKGeometryData<T>& gjk_geometry);

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/gjk_solver_cvx-inl.h"