//
// Created by mech-mind_gw on 3/31/2023.
//

#pragma once

#include "fcl/narrowphase/collision_penetration-inl.h"

namespace fcl {

//==============================================================================
template <typename S>
std::size_t collide(const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
                    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
                    const CollisionRequest<S>& request,
                    CollisionResult<S>& result) {
  // The actual interface
  detail::GJKSolver<S> solver;
  solver.gjk_tolerance = request.binaryCollisionTolerance();
  solver.epa_tolerance = request.distanceTolerance();
  auto penetration_type = request.penetrationMode().penetration_type();
  if (penetration_type ==
          detail::CollisionPenetrationType::DirectedPenetration ||
      penetration_type ==
          detail::CollisionPenetrationType::IncrementalMinimumPenetration) {
    return ::fcl::detail::collisionPenetrationMPR(o1, tf1, o2, tf2, solver,
                                                  request, result);
  } else {
    // No collision penetration or default EPA penetration
    return ::fcl::detail::collide(o1, tf1, o2, tf2, &solver, request, result);
  }
}

//==============================================================================
template <typename S>
std::size_t collide(const CollisionObject<S>* o1, const CollisionObject<S>* o2,
                    const CollisionRequest<S>& request,
                    CollisionResult<S>& result) {
  // Invoke the interface above
  return ::fcl::collide(o1->collisionGeometry().get(), o1->getTransform(),
                        o2->collisionGeometry().get(), o2->getTransform(),
                        request, result);
}

}  // namespace fcl
