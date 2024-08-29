//
// Created by wei on 23-10-16.
//

#pragma once

#include <stack>

#include "fcl/geometry/bvh/BVH_model.h"

namespace fcl {
namespace detail {

template <typename BV>
struct OrientedNodeBVHSolver {
  using S = typename BV::S;
  explicit OrientedNodeBVHSolver(const GJKSolver<S>* gjk_solver_in)
      : gjk_solver(gjk_solver_in) {}
  ~OrientedNodeBVHSolver() = default;

  template <typename Shape>
  void MeshShapeIntersect(const BVHModel<BV>* bvh_1, const Shape& shape_2,
                          const Transform3<S>& tf1, const Transform3<S>& tf2,
                          const CollisionRequest<S>& request,
                          CollisionResult<S>& result) const;
  void MeshIntersect(const BVHModel<BV>* bvh_1, const BVHModel<BV>* bvh_2,
                     const Transform3<S>& tf1, const Transform3<S>& tf2,
                     const CollisionRequest<S>& request,
                     CollisionResult<S>& result) const;

  // Internal state for bvh solver
  const GJKSolver<S>* gjk_solver{nullptr};
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/traversal/collision/bvh_solver-inl.h"