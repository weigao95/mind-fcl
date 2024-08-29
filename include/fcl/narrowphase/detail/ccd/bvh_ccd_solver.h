//
// Created by Wei Gao on 2024/6/14.
//

#pragma once

#include <stack>

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"

namespace fcl {
namespace detail {

template <typename BV>
struct TranslationalDisplacementOrientedNodeBVHSolver {
  using S = typename BV::S;

  template <typename Shape>
  static void RunShapeMesh(const Shape* s1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const BVHModel<BV>* bvh2, const Transform3<S>& tf2,
                           const ContinuousCollisionRequest<S>& request,
                           ContinuousCollisionResult<S>& result);
  template <typename Shape>
  static void RunMeshShape(const BVHModel<BV>* bvh1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const Shape* s2, const Transform3<S>& tf2,
                           const ContinuousCollisionRequest<S>& request,
                           ContinuousCollisionResult<S>& result);

  static void RunMeshPair(const BVHModel<BV>* bvh1, const Transform3<S>& tf1,
                          const TranslationalDisplacement<S>& bvh1_displacement,
                          const BVHModel<BV>* bvh2, const Transform3<S>& tf2,
                          const ContinuousCollisionRequest<S>& request,
                          ContinuousCollisionResult<S>& result);
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/bvh_ccd_solver-inl.h"
