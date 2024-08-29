//
// Created by wei on 24-6-15.
//

#pragma once

#include "fcl/geometry/collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/ccd_request.h"
#include "fcl/narrowphase/detail/ccd/ccd_result.h"

namespace fcl {
namespace detail {

template <typename S>
struct TranslationalCollisionFunctionMatrix {
  /// Typedef
  using CollisionFunc = void (*)(const CollisionGeometry<S>*,
                                 const Transform3<S>&,
                                 const TranslationalDisplacement<S>&,
                                 const CollisionGeometry<S>*,
                                 const Transform3<S>&,
                                 const ContinuousCollisionRequest<S>&,
                                 ContinuousCollisionResult<S>&);

  /// Each item in the collision matrix is a function to handle collision
  /// between objects of type1 and type2
  CollisionFunc collision_matrix[NODE_COUNT][NODE_COUNT];

  /// Construct the interface
  explicit TranslationalCollisionFunctionMatrix();
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/translational_collision_func_matrix-inl.h"