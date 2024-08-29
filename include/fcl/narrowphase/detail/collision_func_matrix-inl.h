/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan */

#ifndef FCL_COLLISION_FUNC_MATRIX_INL_H
#define FCL_COLLISION_FUNC_MATRIX_INL_H

#include "fcl/common/unused.h"
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
#include "fcl/geometry/shape/utility.h"
#include "fcl/narrowphase/collision_object.h"
#include "fcl/narrowphase/detail/collision_func_matrix.h"
#include "fcl/narrowphase/detail/shape_pair_intersect.h"
#include "fcl/narrowphase/detail/traversal/collision/bvh_collision_traversal_node.h"
#include "fcl/narrowphase/detail/traversal/collision/bvh_shape_collision_traversal_node.h"
#include "fcl/narrowphase/detail/traversal/collision/bvh_solver.h"
#include "fcl/narrowphase/detail/traversal/collision/collision_traversal_node_base.h"
#include "fcl/narrowphase/detail/traversal/collision/mesh_collision_traversal_node.h"
#include "fcl/narrowphase/detail/traversal/collision/mesh_shape_collision_traversal_node.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"

// Octree2
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver.h"

// Heightmap
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"

namespace fcl {

namespace detail {

/// Heightmap related method
//==============================================================================
template <typename Shape>
std::size_t ShapeHeightMapCollide(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const GJKSolver<typename Shape::S>* nsolver,
    const CollisionRequest<typename Shape::S>& request,
    CollisionResult<typename Shape::S>& result) {
  using S = typename Shape::S;
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  const auto* obj1 = static_cast<const Shape*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.ShapeHeightMapIntersect(*obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename Shape>
std::size_t HeightMapShapeCollide(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const GJKSolver<typename Shape::S>* nsolver,
    const CollisionRequest<typename Shape::S>& request,
    CollisionResult<typename Shape::S>& result) {
  using S = typename Shape::S;
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  // Perform collision detection
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.HeightMapShapeIntersect(obj1, *obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename S>
std::size_t HeightMapPairCollide(const CollisionGeometry<S>* o1,
                                 const Transform3<S>& tf1,
                                 const CollisionGeometry<S>* o2,
                                 const Transform3<S>& tf2,
                                 const GJKSolver<S>* nsolver,
                                 const CollisionRequest<S>& request,
                                 CollisionResult<S>& result) {
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  // Perform collision detection
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.HeightMapIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename S>
std::size_t HeightMapOctree2Collide(const CollisionGeometry<S>* o1,
                                    const Transform3<S>& tf1,
                                    const CollisionGeometry<S>* o2,
                                    const Transform3<S>& tf2,
                                    const GJKSolver<S>* nsolver,
                                    const CollisionRequest<S>& request,
                                    CollisionResult<S>& result) {
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  // Perform collision detection
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.HeightMapOctreeIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename S>
std::size_t Octree2HeightMapCollide(const CollisionGeometry<S>* o1,
                                    const Transform3<S>& tf1,
                                    const CollisionGeometry<S>* o2,
                                    const Transform3<S>& tf2,
                                    const GJKSolver<S>* nsolver,
                                    const CollisionRequest<S>& request,
                                    CollisionResult<S>& result) {
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  // Perform collision detection
  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.OctreeHeightMapIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename BV>
std::size_t HeightMapBVHCollide(const CollisionGeometry<typename BV::S>* o1,
                                const Transform3<typename BV::S>& tf1,
                                const CollisionGeometry<typename BV::S>* o2,
                                const Transform3<typename BV::S>& tf2,
                                const GJKSolver<typename BV::S>* nsolver,
                                const CollisionRequest<typename BV::S>& request,
                                CollisionResult<typename BV::S>& result) {
  using S = typename BV::S;
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.HeightMapBVHIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename BV>
std::size_t BVHHeightMapCollide(const CollisionGeometry<typename BV::S>* o1,
                                const Transform3<typename BV::S>& tf1,
                                const CollisionGeometry<typename BV::S>* o2,
                                const Transform3<typename BV::S>& tf2,
                                const GJKSolver<typename BV::S>* nsolver,
                                const CollisionRequest<typename BV::S>& request,
                                CollisionResult<typename BV::S>& result) {
  using S = typename BV::S;
  if (request.terminationConditionSatisfied(result)) {
    return result.numContacts();
  }

  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);
  HeightMapCollisionSolver<S> hm_solver(nsolver);
  hm_solver.BVH_HeightMapIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

/// Octree related collision method
//==============================================================================
template <typename Shape>
std::size_t ShapeOcTree2Collide(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const GJKSolver<typename Shape::S>* nsolver,
    const CollisionRequest<typename Shape::S>& request,
    CollisionResult<typename Shape::S>& result) {
  using S = typename Shape::S;
  if (request.terminationConditionSatisfied(result) || o1 == nullptr ||
      o2 == nullptr) {
    return result.numContacts();
  }

  // Obtain the shape
  const auto* obj1 = static_cast<const Shape*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);
  CollisionSolverOctree2<S> octree_solver(nsolver);
  octree_solver.template ShapeOctreeIntersect<Shape>(*obj1, obj2, tf1, tf2,
                                                     request, result);
  return result.numContacts();
}

//==============================================================================
template <typename Shape>
std::size_t OcTree2ShapeCollide(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const GJKSolver<typename Shape::S>* nsolver,
    const CollisionRequest<typename Shape::S>& request,
    CollisionResult<typename Shape::S>& result) {
  using S = typename Shape::S;
  if (request.terminationConditionSatisfied(result) || o1 == nullptr ||
      o2 == nullptr) {
    return result.numContacts();
  }

  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);
  CollisionSolverOctree2<S> octree_solver(nsolver);
  octree_solver.template OctreeShapeIntersect<Shape>(obj1, *obj2, tf1, tf2,
                                                     request, result);
  return result.numContacts();
}

//==============================================================================
template <typename S>
std::size_t OcTree2Collide(const CollisionGeometry<S>* o1,
                           const Transform3<S>& tf1,
                           const CollisionGeometry<S>* o2,
                           const Transform3<S>& tf2,
                           const GJKSolver<S>* nsolver,
                           const CollisionRequest<S>& request,
                           CollisionResult<S>& result) {
  if (request.terminationConditionSatisfied(result) || o1 == nullptr ||
      o2 == nullptr)
    return result.numContacts();

  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);
  CollisionSolverOctree2<S> octree_solver(nsolver);
  octree_solver.OctreePairIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename BV>
std::size_t OcTree2BVHCollide(const CollisionGeometry<typename BV::S>* o1,
                              const Transform3<typename BV::S>& tf1,
                              const CollisionGeometry<typename BV::S>* o2,
                              const Transform3<typename BV::S>& tf2,
                              const GJKSolver<typename BV::S>* nsolver,
                              const CollisionRequest<typename BV::S>& request,
                              CollisionResult<typename BV::S>& result) {
  using S = typename BV::S;
  if (request.terminationConditionSatisfied(result))
    return result.numContacts();

  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);
  CollisionSolverOctree2<S> octree_solver(nsolver);
  octree_solver.template OctreeBVHIntersect<BV>(obj1, obj2, tf1, tf2, request,
                                                result);
  return result.numContacts();
}

//==============================================================================
template <typename BV>
std::size_t BVHOcTree2Collide(const CollisionGeometry<typename BV::S>* o1,
                              const Transform3<typename BV::S>& tf1,
                              const CollisionGeometry<typename BV::S>* o2,
                              const Transform3<typename BV::S>& tf2,
                              const GJKSolver<typename BV::S>* nsolver,
                              const CollisionRequest<typename BV::S>& request,
                              CollisionResult<typename BV::S>& result) {
  using S = typename BV::S;
  if (request.terminationConditionSatisfied(result))
    return result.numContacts();

  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);
  CollisionSolverOctree2<S> octree_solver(nsolver);
  octree_solver.template BVH_OctreeIntersect<BV>(obj1, obj2, tf1, tf2, request,
                                                 result);
  return result.numContacts();
}

/// Shape and BVH
//==============================================================================
template <typename Shape1, typename Shape2>
std::size_t ShapeShapeCollide(
    const CollisionGeometry<typename Shape1::S>* o1,
    const Transform3<typename Shape1::S>& tf1,
    const CollisionGeometry<typename Shape1::S>* o2,
    const Transform3<typename Shape1::S>& tf2,
    const GJKSolver<typename Shape1::S>* nsolver,
    const CollisionRequest<typename Shape1::S>& request,
    CollisionResult<typename Shape1::S>& result) {
  if (request.terminationConditionSatisfied(result))
    return result.numContacts();

  // Obtain the objects
  using S = typename Shape1::S;
  const auto* obj1 = static_cast<const Shape1*>(o1);
  const auto* obj2 = static_cast<const Shape2*>(o2);
  ShapePairIntersectSolver<S> shape_solver(nsolver);
  shape_solver.template ShapeIntersect<Shape1, Shape2>(obj1, tf1, obj2, tf2,
                                                       request, result);
  return result.numContacts();
}

//==============================================================================
template <typename BV, typename Shape>
struct BVHShapeCollider {
  using S = typename BV::S;

  static std::size_t collide(const CollisionGeometry<S>* o1,
                             const Transform3<S>& tf1,
                             const CollisionGeometry<S>* o2,
                             const Transform3<S>& tf2,
                             const GJKSolver<S>* nsolver,
                             const CollisionRequest<S>& request,
                             CollisionResult<S>& result) {
    if (request.terminationConditionSatisfied(result))
      return result.numContacts();

    // Copy the data to avoid mutation
    const BVHModel<BV>* obj1 = static_cast<const BVHModel<BV>*>(o1);
    std::unique_ptr<BVHModel<BV>> obj1_tmp(new BVHModel<BV>(*obj1));
    Transform3<S> tf1_tmp = tf1;
    const Shape* obj2 = static_cast<const Shape*>(o2);

    MeshShapeCollisionTraversalNode<BV, Shape> node;
    initialize(node, *obj1_tmp, tf1_tmp, *obj2, tf2, nsolver, request, result);
    fcl::detail::collide(&node);
    return result.numContacts();
  }
};

//==============================================================================
template <typename BV, typename Shape>
std::size_t orientedBVHShapeCollide(
    const CollisionGeometry<typename BV::S>* o1,
    const Transform3<typename BV::S>& tf1,
    const CollisionGeometry<typename BV::S>* o2,
    const Transform3<typename BV::S>& tf2,
    const GJKSolver<typename BV::S>* nsolver,
    const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result) {
  if (request.terminationConditionSatisfied(result))
    return result.numContacts();

  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);
  OrientedNodeBVHSolver<BV> bvh_solver(nsolver);
  bvh_solver.template MeshShapeIntersect<Shape>(obj1, *obj2, tf1, tf2, request,
                                                result);
  return result.numContacts();
}

//==============================================================================
template <typename Shape>
struct BVHShapeCollider<OBB<typename Shape::S>, Shape> {
  using S = typename Shape::S;

  static std::size_t collide(const CollisionGeometry<S>* o1,
                             const Transform3<S>& tf1,
                             const CollisionGeometry<S>* o2,
                             const Transform3<S>& tf2,
                             const GJKSolver<S>* nsolver,
                             const CollisionRequest<S>& request,
                             CollisionResult<S>& result) {
    return detail::orientedBVHShapeCollide<OBB<S>, Shape>(
        o1, tf1, o2, tf2, nsolver, request, result);
  }
};

//==============================================================================
template <typename Shape>
struct BVHShapeCollider<RSS<typename Shape::S>, Shape> {
  using S = typename Shape::S;

  static std::size_t collide(const CollisionGeometry<S>* o1,
                             const Transform3<S>& tf1,
                             const CollisionGeometry<S>* o2,
                             const Transform3<S>& tf2,
                             const GJKSolver<S>* nsolver,
                             const CollisionRequest<S>& request,
                             CollisionResult<S>& result) {
    return detail::orientedBVHShapeCollide<RSS<S>, Shape>(
        o1, tf1, o2, tf2, nsolver, request, result);
  }
};

//==============================================================================
template <typename Shape>
struct BVHShapeCollider<kIOS<typename Shape::S>, Shape> {
  using S = typename Shape::S;

  static std::size_t collide(const CollisionGeometry<S>* o1,
                             const Transform3<S>& tf1,
                             const CollisionGeometry<S>* o2,
                             const Transform3<S>& tf2,
                             const GJKSolver<S>* nsolver,
                             const CollisionRequest<S>& request,
                             CollisionResult<S>& result) {
    return detail::orientedBVHShapeCollide<kIOS<S>, Shape>(
        o1, tf1, o2, tf2, nsolver, request, result);
  }
};

//==============================================================================
template <typename Shape>
struct BVHShapeCollider<OBBRSS<typename Shape::S>, Shape> {
  using S = typename Shape::S;

  static std::size_t collide(const CollisionGeometry<S>* o1,
                             const Transform3<S>& tf1,
                             const CollisionGeometry<S>* o2,
                             const Transform3<S>& tf2,
                             const GJKSolver<S>* nsolver,
                             const CollisionRequest<S>& request,
                             CollisionResult<S>& result) {
    return detail::orientedBVHShapeCollide<OBBRSS<S>, Shape>(
        o1, tf1, o2, tf2, nsolver, request, result);
  }
};

//==============================================================================
template <typename S, typename BV>
struct BVHCollideImpl {
  static std::size_t run(const CollisionGeometry<S>* o1,
                         const Transform3<S>& tf1,
                         const CollisionGeometry<S>* o2,
                         const Transform3<S>& tf2, const GJKSolver<S>* nsolver,
                         const CollisionRequest<S>& request,
                         CollisionResult<S>& result) {
    if (request.terminationConditionSatisfied(result))
      return result.numContacts();

    // Copy the data to avoid mutation
    assert(nsolver != nullptr);
    const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
    const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);
    std::unique_ptr<BVHModel<BV>> obj1_tmp(new BVHModel<BV>(*obj1));
    Transform3<S> tf1_tmp = tf1;
    std::unique_ptr<BVHModel<BV>> obj2_tmp(new BVHModel<BV>(*obj2));
    Transform3<S> tf2_tmp = tf2;

    MeshCollisionTraversalNode<BV> node;
    initialize(node, *obj1_tmp, tf1_tmp, *obj2_tmp, tf2_tmp, nsolver, request,
               result);
    collide(&node);
    return result.numContacts();
  }
};

//==============================================================================
template <typename BV>
std::size_t orientedMeshCollide(const CollisionGeometry<typename BV::S>* o1,
                                const Transform3<typename BV::S>& tf1,
                                const CollisionGeometry<typename BV::S>* o2,
                                const Transform3<typename BV::S>& tf2,
                                const GJKSolver<typename BV::S>* nsolver,
                                const CollisionRequest<typename BV::S>& request,
                                CollisionResult<typename BV::S>& result) {
  if (request.terminationConditionSatisfied(result))
    return result.numContacts();

  assert(nsolver != nullptr);
  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);

  OrientedNodeBVHSolver<BV> bvh_solver(nsolver);
  bvh_solver.MeshIntersect(obj1, obj2, tf1, tf2, request, result);
  return result.numContacts();
}

//==============================================================================
template <typename S>
struct BVHCollideImpl<S, OBB<S>> {
  static std::size_t run(const CollisionGeometry<S>* o1,
                         const Transform3<S>& tf1,
                         const CollisionGeometry<S>* o2,
                         const Transform3<S>& tf2, const GJKSolver<S>* nsolver,
                         const CollisionRequest<S>& request,
                         CollisionResult<S>& result) {
    return detail::orientedMeshCollide<OBB<S>>(o1, tf1, o2, tf2, nsolver,
                                               request, result);
  }
};

//==============================================================================
template <typename S>
struct BVHCollideImpl<S, OBBRSS<S>> {
  static std::size_t run(const CollisionGeometry<S>* o1,
                         const Transform3<S>& tf1,
                         const CollisionGeometry<S>* o2,
                         const Transform3<S>& tf2, const GJKSolver<S>* nsolver,
                         const CollisionRequest<S>& request,
                         CollisionResult<S>& result) {
    return detail::orientedMeshCollide<OBBRSS<S>>(o1, tf1, o2, tf2, nsolver,
                                                  request, result);
  }
};

//==============================================================================
template <typename S>
struct BVHCollideImpl<S, kIOS<S>> {
  static std::size_t run(const CollisionGeometry<S>* o1,
                         const Transform3<S>& tf1,
                         const CollisionGeometry<S>* o2,
                         const Transform3<S>& tf2, const GJKSolver<S>* nsolver,
                         const CollisionRequest<S>& request,
                         CollisionResult<S>& result) {
    return detail::orientedMeshCollide<kIOS<S>>(o1, tf1, o2, tf2, nsolver,
                                                request, result);
  }
};

//==============================================================================
template <typename BV>
std::size_t BVHCollide(const CollisionGeometry<typename BV::S>* o1,
                       const Transform3<typename BV::S>& tf1,
                       const CollisionGeometry<typename BV::S>* o2,
                       const Transform3<typename BV::S>& tf2,
                       const GJKSolver<typename BV::S>* nsolver,
                       const CollisionRequest<typename BV::S>& request,
                       CollisionResult<typename BV::S>& result) {
  // return BVHCollide<BV>(o1, tf1, o2, tf2, request, result);
  return BVHCollideImpl<typename BV::S, BV>::run(o1, tf1, o2, tf2, nsolver,
                                                 request, result);
}

//==============================================================================
template <typename S_>
CollisionFunctionMatrix<S_>::CollisionFunctionMatrix() {
  for (int i = 0; i < NODE_COUNT; ++i) {
    for (int j = 0; j < NODE_COUNT; ++j) collision_matrix[i][j] = nullptr;
  }

  /// Shape vs shape
  // clang-format off
  collision_matrix[GEOM_BOX][GEOM_BOX] = &ShapeShapeCollide<Box<S>, Box<S>>;
  collision_matrix[GEOM_BOX][GEOM_SPHERE] =&ShapeShapeCollide<Box<S>, Sphere<S>>;
  collision_matrix[GEOM_BOX][GEOM_ELLIPSOID] =&ShapeShapeCollide<Box<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_BOX][GEOM_CAPSULE] =&ShapeShapeCollide<Box<S>, Capsule<S>>;
  collision_matrix[GEOM_BOX][GEOM_CONE] = &ShapeShapeCollide<Box<S>, Cone<S>>;
  collision_matrix[GEOM_BOX][GEOM_CYLINDER] =&ShapeShapeCollide<Box<S>, Cylinder<S>>;
  collision_matrix[GEOM_BOX][GEOM_CONVEX] =&ShapeShapeCollide<Box<S>, Convex<S>>;
  collision_matrix[GEOM_BOX][GEOM_PLANE] = &ShapeShapeCollide<Box<S>, Plane<S>>;
  collision_matrix[GEOM_BOX][GEOM_HALFSPACE] =&ShapeShapeCollide<Box<S>, Halfspace<S>>;

  collision_matrix[GEOM_SPHERE][GEOM_BOX] =&ShapeShapeCollide<Sphere<S>, Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_SPHERE] =&ShapeShapeCollide<Sphere<S>, Sphere<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_ELLIPSOID] =&ShapeShapeCollide<Sphere<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CAPSULE] =&ShapeShapeCollide<Sphere<S>, Capsule<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CONE] =&ShapeShapeCollide<Sphere<S>, Cone<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CYLINDER] =&ShapeShapeCollide<Sphere<S>, Cylinder<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CONVEX] =&ShapeShapeCollide<Sphere<S>, Convex<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_PLANE] =&ShapeShapeCollide<Sphere<S>, Plane<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_HALFSPACE] =&ShapeShapeCollide<Sphere<S>, Halfspace<S>>;

  collision_matrix[GEOM_ELLIPSOID][GEOM_BOX] =&ShapeShapeCollide<Ellipsoid<S>, Box<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_SPHERE] =&ShapeShapeCollide<Ellipsoid<S>, Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_ELLIPSOID] =&ShapeShapeCollide<Ellipsoid<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CAPSULE] =&ShapeShapeCollide<Ellipsoid<S>, Capsule<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CONE] =&ShapeShapeCollide<Ellipsoid<S>, Cone<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CYLINDER] =&ShapeShapeCollide<Ellipsoid<S>, Cylinder<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CONVEX] =&ShapeShapeCollide<Ellipsoid<S>, Convex<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_PLANE] =&ShapeShapeCollide<Ellipsoid<S>, Plane<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_HALFSPACE] =&ShapeShapeCollide<Ellipsoid<S>, Halfspace<S>>;

  collision_matrix[GEOM_CAPSULE][GEOM_BOX] =&ShapeShapeCollide<Capsule<S>, Box<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_SPHERE] =&ShapeShapeCollide<Capsule<S>, Sphere<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_ELLIPSOID] =&ShapeShapeCollide<Capsule<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CAPSULE] =&ShapeShapeCollide<Capsule<S>, Capsule<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CONE] =&ShapeShapeCollide<Capsule<S>, Cone<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CYLINDER] =&ShapeShapeCollide<Capsule<S>, Cylinder<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CONVEX] =&ShapeShapeCollide<Capsule<S>, Convex<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_PLANE] =&ShapeShapeCollide<Capsule<S>, Plane<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_HALFSPACE] =&ShapeShapeCollide<Capsule<S>, Halfspace<S>>;

  collision_matrix[GEOM_CONE][GEOM_BOX] = &ShapeShapeCollide<Cone<S>, Box<S>>;
  collision_matrix[GEOM_CONE][GEOM_SPHERE] =&ShapeShapeCollide<Cone<S>, Sphere<S>>;
  collision_matrix[GEOM_CONE][GEOM_ELLIPSOID] =&ShapeShapeCollide<Cone<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CONE][GEOM_CAPSULE] =&ShapeShapeCollide<Cone<S>, Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_CONE] = &ShapeShapeCollide<Cone<S>, Cone<S>>;
  collision_matrix[GEOM_CONE][GEOM_CYLINDER] =&ShapeShapeCollide<Cone<S>, Cylinder<S>>;
  collision_matrix[GEOM_CONE][GEOM_CONVEX] =&ShapeShapeCollide<Cone<S>, Convex<S>>;
  collision_matrix[GEOM_CONE][GEOM_PLANE] =&ShapeShapeCollide<Cone<S>, Plane<S>>;
  collision_matrix[GEOM_CONE][GEOM_HALFSPACE] =&ShapeShapeCollide<Cone<S>, Halfspace<S>>;

  collision_matrix[GEOM_CYLINDER][GEOM_BOX] =&ShapeShapeCollide<Cylinder<S>, Box<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_SPHERE] =&ShapeShapeCollide<Cylinder<S>, Sphere<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_ELLIPSOID] =&ShapeShapeCollide<Cylinder<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CAPSULE] =&ShapeShapeCollide<Cylinder<S>, Capsule<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CONE] =&ShapeShapeCollide<Cylinder<S>, Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CYLINDER] =&ShapeShapeCollide<Cylinder<S>, Cylinder<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CONVEX] =&ShapeShapeCollide<Cylinder<S>, Convex<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_PLANE] =&ShapeShapeCollide<Cylinder<S>, Plane<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_HALFSPACE] =&ShapeShapeCollide<Cylinder<S>, Halfspace<S>>;

  collision_matrix[GEOM_CONVEX][GEOM_BOX] =&ShapeShapeCollide<Convex<S>, Box<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_SPHERE] =&ShapeShapeCollide<Convex<S>, Sphere<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_ELLIPSOID] =&ShapeShapeCollide<Convex<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CAPSULE] =&ShapeShapeCollide<Convex<S>, Capsule<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CONE] =&ShapeShapeCollide<Convex<S>, Cone<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CYLINDER] =&ShapeShapeCollide<Convex<S>, Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CONVEX] =&ShapeShapeCollide<Convex<S>, Convex<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_PLANE] =&ShapeShapeCollide<Convex<S>, Plane<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_HALFSPACE] =&ShapeShapeCollide<Convex<S>, Halfspace<S>>;

  collision_matrix[GEOM_PLANE][GEOM_BOX] = &ShapeShapeCollide<Plane<S>, Box<S>>;
  collision_matrix[GEOM_PLANE][GEOM_SPHERE] =&ShapeShapeCollide<Plane<S>, Sphere<S>>;
  collision_matrix[GEOM_PLANE][GEOM_ELLIPSOID] =&ShapeShapeCollide<Plane<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CAPSULE] =&ShapeShapeCollide<Plane<S>, Capsule<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CONE] =&ShapeShapeCollide<Plane<S>, Cone<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CYLINDER] =&ShapeShapeCollide<Plane<S>, Cylinder<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CONVEX] =&ShapeShapeCollide<Plane<S>, Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_PLANE] =&ShapeShapeCollide<Plane<S>, Plane<S>>;
  collision_matrix[GEOM_PLANE][GEOM_HALFSPACE] =&ShapeShapeCollide<Plane<S>, Halfspace<S>>;

  collision_matrix[GEOM_HALFSPACE][GEOM_BOX] =&ShapeShapeCollide<Halfspace<S>, Box<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_SPHERE] =&ShapeShapeCollide<Halfspace<S>, Sphere<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_ELLIPSOID] =&ShapeShapeCollide<Halfspace<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CAPSULE] =&ShapeShapeCollide<Halfspace<S>, Capsule<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CONE] =&ShapeShapeCollide<Halfspace<S>, Cone<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CYLINDER] =&ShapeShapeCollide<Halfspace<S>, Cylinder<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CONVEX] =&ShapeShapeCollide<Halfspace<S>, Convex<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_PLANE] =&ShapeShapeCollide<Halfspace<S>, Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_HALFSPACE] =&ShapeShapeCollide<Halfspace<S>, Halfspace<S>>;

  /// BVH
  collision_matrix[BV_AABB][GEOM_BOX] =&BVHShapeCollider<AABB<S>, Box<S>>::collide;
  collision_matrix[BV_AABB][GEOM_SPHERE] =&BVHShapeCollider<AABB<S>, Sphere<S>>::collide;
  collision_matrix[BV_AABB][GEOM_ELLIPSOID] =&BVHShapeCollider<AABB<S>, Ellipsoid<S>>::collide;
  collision_matrix[BV_AABB][GEOM_CAPSULE] =&BVHShapeCollider<AABB<S>, Capsule<S>>::collide;
  collision_matrix[BV_AABB][GEOM_CONE] =&BVHShapeCollider<AABB<S>, Cone<S>>::collide;
  collision_matrix[BV_AABB][GEOM_CYLINDER] =&BVHShapeCollider<AABB<S>, Cylinder<S>>::collide;
  collision_matrix[BV_AABB][GEOM_CONVEX] =&BVHShapeCollider<AABB<S>, Convex<S>>::collide;
  collision_matrix[BV_AABB][GEOM_PLANE] =&BVHShapeCollider<AABB<S>, Plane<S>>::collide;
  collision_matrix[BV_AABB][GEOM_HALFSPACE] =&BVHShapeCollider<AABB<S>, Halfspace<S>>::collide;

  collision_matrix[BV_OBB][GEOM_BOX] =&BVHShapeCollider<OBB<S>, Box<S>>::collide;
  collision_matrix[BV_OBB][GEOM_SPHERE] =&BVHShapeCollider<OBB<S>, Sphere<S>>::collide;
  collision_matrix[BV_OBB][GEOM_ELLIPSOID] =&BVHShapeCollider<OBB<S>, Ellipsoid<S>>::collide;
  collision_matrix[BV_OBB][GEOM_CAPSULE] =&BVHShapeCollider<OBB<S>, Capsule<S>>::collide;
  collision_matrix[BV_OBB][GEOM_CONE] =&BVHShapeCollider<OBB<S>, Cone<S>>::collide;
  collision_matrix[BV_OBB][GEOM_CYLINDER] =&BVHShapeCollider<OBB<S>, Cylinder<S>>::collide;
  collision_matrix[BV_OBB][GEOM_CONVEX] =&BVHShapeCollider<OBB<S>, Convex<S>>::collide;
  collision_matrix[BV_OBB][GEOM_PLANE] =&BVHShapeCollider<OBB<S>, Plane<S>>::collide;
  collision_matrix[BV_OBB][GEOM_HALFSPACE] =&BVHShapeCollider<OBB<S>, Halfspace<S>>::collide;

  collision_matrix[BV_RSS][GEOM_BOX] =&BVHShapeCollider<RSS<S>, Box<S>>::collide;
  collision_matrix[BV_RSS][GEOM_SPHERE] =&BVHShapeCollider<RSS<S>, Sphere<S>>::collide;
  collision_matrix[BV_RSS][GEOM_ELLIPSOID] =&BVHShapeCollider<RSS<S>, Ellipsoid<S>>::collide;
  collision_matrix[BV_RSS][GEOM_CAPSULE] =&BVHShapeCollider<RSS<S>, Capsule<S>>::collide;
  collision_matrix[BV_RSS][GEOM_CONE] =&BVHShapeCollider<RSS<S>, Cone<S>>::collide;
  collision_matrix[BV_RSS][GEOM_CYLINDER] =&BVHShapeCollider<RSS<S>, Cylinder<S>>::collide;
  collision_matrix[BV_RSS][GEOM_CONVEX] =&BVHShapeCollider<RSS<S>, Convex<S>>::collide;
  collision_matrix[BV_RSS][GEOM_PLANE] =&BVHShapeCollider<RSS<S>, Plane<S>>::collide;
  collision_matrix[BV_RSS][GEOM_HALFSPACE] =&BVHShapeCollider<RSS<S>, Halfspace<S>>::collide;

  collision_matrix[BV_KDOP16][GEOM_BOX] =&BVHShapeCollider<KDOP<S, 16>, Box<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_SPHERE] =&BVHShapeCollider<KDOP<S, 16>, Sphere<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_ELLIPSOID] =&BVHShapeCollider<KDOP<S, 16>, Ellipsoid<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_CAPSULE] =&BVHShapeCollider<KDOP<S, 16>, Capsule<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_CONE] =&BVHShapeCollider<KDOP<S, 16>, Cone<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_CYLINDER] =&BVHShapeCollider<KDOP<S, 16>, Cylinder<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_CONVEX] =&BVHShapeCollider<KDOP<S, 16>, Convex<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_PLANE] =&BVHShapeCollider<KDOP<S, 16>, Plane<S>>::collide;
  collision_matrix[BV_KDOP16][GEOM_HALFSPACE] =&BVHShapeCollider<KDOP<S, 16>, Halfspace<S>>::collide;

  collision_matrix[BV_KDOP18][GEOM_BOX] =&BVHShapeCollider<KDOP<S, 18>, Box<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_SPHERE] =&BVHShapeCollider<KDOP<S, 18>, Sphere<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_ELLIPSOID] =&BVHShapeCollider<KDOP<S, 18>, Ellipsoid<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_CAPSULE] =&BVHShapeCollider<KDOP<S, 18>, Capsule<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_CONE] =&BVHShapeCollider<KDOP<S, 18>, Cone<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_CYLINDER] =&BVHShapeCollider<KDOP<S, 18>, Cylinder<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_CONVEX] =&BVHShapeCollider<KDOP<S, 18>, Convex<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_PLANE] =&BVHShapeCollider<KDOP<S, 18>, Plane<S>>::collide;
  collision_matrix[BV_KDOP18][GEOM_HALFSPACE] =&BVHShapeCollider<KDOP<S, 18>, Halfspace<S>>::collide;

  collision_matrix[BV_KDOP24][GEOM_BOX] =&BVHShapeCollider<KDOP<S, 24>, Box<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_SPHERE] =&BVHShapeCollider<KDOP<S, 24>, Sphere<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_ELLIPSOID] =&BVHShapeCollider<KDOP<S, 24>, Ellipsoid<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_CAPSULE] =&BVHShapeCollider<KDOP<S, 24>, Capsule<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_CONE] =&BVHShapeCollider<KDOP<S, 24>, Cone<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_CYLINDER] =&BVHShapeCollider<KDOP<S, 24>, Cylinder<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_CONVEX] =&BVHShapeCollider<KDOP<S, 24>, Convex<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_PLANE] =&BVHShapeCollider<KDOP<S, 24>, Plane<S>>::collide;
  collision_matrix[BV_KDOP24][GEOM_HALFSPACE] =&BVHShapeCollider<KDOP<S, 24>, Halfspace<S>>::collide;

  collision_matrix[BV_kIOS][GEOM_BOX] =&BVHShapeCollider<kIOS<S>, Box<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_SPHERE] =&BVHShapeCollider<kIOS<S>, Sphere<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_ELLIPSOID] =&BVHShapeCollider<kIOS<S>, Ellipsoid<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_CAPSULE] =&BVHShapeCollider<kIOS<S>, Capsule<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_CONE] =&BVHShapeCollider<kIOS<S>, Cone<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_CYLINDER] =&BVHShapeCollider<kIOS<S>, Cylinder<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_CONVEX] =&BVHShapeCollider<kIOS<S>, Convex<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_PLANE] =&BVHShapeCollider<kIOS<S>, Plane<S>>::collide;
  collision_matrix[BV_kIOS][GEOM_HALFSPACE] =&BVHShapeCollider<kIOS<S>, Halfspace<S>>::collide;

  collision_matrix[BV_OBBRSS][GEOM_BOX] =&BVHShapeCollider<OBBRSS<S>, Box<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_SPHERE] =&BVHShapeCollider<OBBRSS<S>, Sphere<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_ELLIPSOID] =&BVHShapeCollider<OBBRSS<S>, Ellipsoid<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_CAPSULE] =&BVHShapeCollider<OBBRSS<S>, Capsule<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_CONE] =&BVHShapeCollider<OBBRSS<S>, Cone<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_CYLINDER] =&BVHShapeCollider<OBBRSS<S>, Cylinder<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_CONVEX] =&BVHShapeCollider<OBBRSS<S>, Convex<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_PLANE] =&BVHShapeCollider<OBBRSS<S>, Plane<S>>::collide;
  collision_matrix[BV_OBBRSS][GEOM_HALFSPACE] =&BVHShapeCollider<OBBRSS<S>, Halfspace<S>>::collide;

  collision_matrix[BV_AABB][BV_AABB] = &BVHCollide<AABB<S>>;
  collision_matrix[BV_OBB][BV_OBB] = &BVHCollide<OBB<S>>;
  collision_matrix[BV_RSS][BV_RSS] = &BVHCollide<RSS<S>>;
  collision_matrix[BV_KDOP16][BV_KDOP16] = &BVHCollide<KDOP<S, 16>>;
  collision_matrix[BV_KDOP18][BV_KDOP18] = &BVHCollide<KDOP<S, 18>>;
  collision_matrix[BV_KDOP24][BV_KDOP24] = &BVHCollide<KDOP<S, 24>>;
  collision_matrix[BV_kIOS][BV_kIOS] = &BVHCollide<kIOS<S>>;
  collision_matrix[BV_OBBRSS][BV_OBBRSS] = &BVHCollide<OBBRSS<S>>;

  /// Octree2
  collision_matrix[GEOM_OCTREE2][GEOM_BOX] = &OcTree2ShapeCollide<Box<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_SPHERE] = &OcTree2ShapeCollide<Sphere<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_ELLIPSOID] = &OcTree2ShapeCollide<Ellipsoid<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CAPSULE] = &OcTree2ShapeCollide<Capsule<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CONE] = &OcTree2ShapeCollide<Cone<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CYLINDER] = &OcTree2ShapeCollide<Cylinder<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CONVEX] = &OcTree2ShapeCollide<Convex<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_PLANE] = &OcTree2ShapeCollide<Plane<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_HALFSPACE] = &OcTree2ShapeCollide<Halfspace<S>>;

  collision_matrix[GEOM_BOX][GEOM_OCTREE2] = &ShapeOcTree2Collide<Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_OCTREE2] = &ShapeOcTree2Collide<Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_OCTREE2] = &ShapeOcTree2Collide<Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_OCTREE2] = &ShapeOcTree2Collide<Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_OCTREE2] = &ShapeOcTree2Collide<Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_OCTREE2] = &ShapeOcTree2Collide<Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_OCTREE2] = &ShapeOcTree2Collide<Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_OCTREE2] = &ShapeOcTree2Collide<Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_OCTREE2] = &ShapeOcTree2Collide<Halfspace<S>>;

  collision_matrix[GEOM_OCTREE2][GEOM_OCTREE2] = &OcTree2Collide<S>;

  collision_matrix[GEOM_OCTREE2][BV_AABB] = &OcTree2BVHCollide<AABB<S>>;
  collision_matrix[GEOM_OCTREE2][BV_OBB] = &OcTree2BVHCollide<OBB<S>>;
  collision_matrix[GEOM_OCTREE2][BV_RSS] = &OcTree2BVHCollide<RSS<S>>;
  collision_matrix[GEOM_OCTREE2][BV_OBBRSS] = &OcTree2BVHCollide<OBBRSS<S>>;
  collision_matrix[GEOM_OCTREE2][BV_kIOS] = &OcTree2BVHCollide<kIOS<S>>;
  collision_matrix[GEOM_OCTREE2][BV_KDOP16] = &OcTree2BVHCollide<KDOP<S, 16>>;
  collision_matrix[GEOM_OCTREE2][BV_KDOP18] = &OcTree2BVHCollide<KDOP<S, 18>>;
  collision_matrix[GEOM_OCTREE2][BV_KDOP24] = &OcTree2BVHCollide<KDOP<S, 24>>;

  collision_matrix[BV_AABB][GEOM_OCTREE2] = &BVHOcTree2Collide<AABB<S>>;
  collision_matrix[BV_OBB][GEOM_OCTREE2] = &BVHOcTree2Collide<OBB<S>>;
  collision_matrix[BV_RSS][GEOM_OCTREE2] = &BVHOcTree2Collide<RSS<S>>;
  collision_matrix[BV_OBBRSS][GEOM_OCTREE2] = &BVHOcTree2Collide<OBBRSS<S>>;
  collision_matrix[BV_kIOS][GEOM_OCTREE2] = &BVHOcTree2Collide<kIOS<S>>;
  collision_matrix[BV_KDOP16][GEOM_OCTREE2] = &BVHOcTree2Collide<KDOP<S, 16>>;
  collision_matrix[BV_KDOP18][GEOM_OCTREE2] = &BVHOcTree2Collide<KDOP<S, 18>>;
  collision_matrix[BV_KDOP24][GEOM_OCTREE2] = &BVHOcTree2Collide<KDOP<S, 24>>;

  /// Heightmap
  collision_matrix[GEOM_HEIGHTMAP][GEOM_BOX] = &HeightMapShapeCollide<Box<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_SPHERE] = &HeightMapShapeCollide<Sphere<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_ELLIPSOID] = &HeightMapShapeCollide<Ellipsoid<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CAPSULE] = &HeightMapShapeCollide<Capsule<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CONE] = &HeightMapShapeCollide<Cone<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CYLINDER] = &HeightMapShapeCollide<Cylinder<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CONVEX] = &HeightMapShapeCollide<Convex<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_PLANE] = &HeightMapShapeCollide<Plane<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_HALFSPACE] = &HeightMapShapeCollide<Halfspace<S>>;

  collision_matrix[GEOM_BOX][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_HEIGHTMAP] = &ShapeHeightMapCollide<Halfspace<S>>;

  collision_matrix[GEOM_HEIGHTMAP][GEOM_HEIGHTMAP] = &HeightMapPairCollide<S>;

  collision_matrix[GEOM_HEIGHTMAP][GEOM_OCTREE2] = &HeightMapOctree2Collide<S>;
  collision_matrix[GEOM_OCTREE2][GEOM_HEIGHTMAP] = &Octree2HeightMapCollide<S>;

  collision_matrix[GEOM_HEIGHTMAP][BV_AABB] = &HeightMapBVHCollide<AABB<S>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_OBB] = &HeightMapBVHCollide<OBB<S>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_RSS] = &HeightMapBVHCollide<RSS<S>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_OBBRSS] = &HeightMapBVHCollide<OBBRSS<S>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_kIOS] = &HeightMapBVHCollide<kIOS<S>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_KDOP16] = &HeightMapBVHCollide<KDOP<S, 16>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_KDOP18] = &HeightMapBVHCollide<KDOP<S, 18>>;
  collision_matrix[GEOM_HEIGHTMAP][BV_KDOP24] = &HeightMapBVHCollide<KDOP<S, 24>>;

  collision_matrix[BV_AABB][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<AABB<S>>;
  collision_matrix[BV_OBB][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<OBB<S>>;
  collision_matrix[BV_RSS][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<RSS<S>>;
  collision_matrix[BV_OBBRSS][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<OBBRSS<S>>;
  collision_matrix[BV_kIOS][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<kIOS<S>>;
  collision_matrix[BV_KDOP16][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<KDOP<S, 16>>;
  collision_matrix[BV_KDOP18][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<KDOP<S, 18>>;
  collision_matrix[BV_KDOP24][GEOM_HEIGHTMAP] = &BVHHeightMapCollide<KDOP<S, 24>>;
  // clang-format on
}

}  // namespace detail
}  // namespace fcl

#endif
