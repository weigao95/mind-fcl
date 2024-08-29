//
// Created by wei on 24-6-15.
//

#pragma once

#include "fcl/geometry/collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/ccd_request.h"
#include "fcl/narrowphase/detail/ccd/ccd_result.h"

namespace fcl {

template <typename S>
void translational_ccd(const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
                       const TranslationalDisplacement<S>& o1_displacement,
                       const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
                       const ContinuousCollisionRequest<S>& request,
                       ContinuousCollisionResult<S>& result);

} // namespace fcl

#include "fcl/narrowphase/continuous_collision-inl.h"