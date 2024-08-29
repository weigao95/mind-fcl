//
// Created by wei on 24-6-15.
//

#pragma once

#include "fcl/geometry/collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"

namespace fcl {
namespace detail {

template <typename Shape>
void initializeShapeFixedOrientationBoxTranslationalCCD(
    const Shape& shape1, const Transform3<typename Shape::S>& tf_shape1,
    const TranslationalDisplacement<typename Shape::S>& displacement_shape1,
    const Transform3<typename Shape::S>& tf_another_2,
    // output
    FixedOrientationBoxPairTranslationalCCD<typename Shape::S>& disjoint,
    AABB<typename Shape::S>& shape_local_AABB);

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/ccd_solver_utility-inl.h"