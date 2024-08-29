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

#ifndef FCL_TRAVERSAL_MESHSHAPECOLLISIONTRAVERSALNODE_H
#define FCL_TRAVERSAL_MESHSHAPECOLLISIONTRAVERSALNODE_H

#include "fcl/geometry/shape/utility.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/traversal/collision/bvh_shape_collision_traversal_node.h"

namespace fcl {

namespace detail {

/// @brief Traversal node for collision between mesh and shape
template <typename BV, typename Shape>
class MeshShapeCollisionTraversalNode
    : public BVHShapeCollisionTraversalNode<BV, Shape> {
 public:
  using S = typename BV::S;
  MeshShapeCollisionTraversalNode();
  ~MeshShapeCollisionTraversalNode() override = default;

  /// Implement the interface
  void leafTesting(int b1, int b2) const override;
  bool canStop() const override;

  S cost_density;
  const GJKSolver<S>* nsolver;
};

/// @brief Initialize traversal node for collision between one mesh and one
/// shape, given current object transform
template <typename BV, typename Shape>
bool initialize(MeshShapeCollisionTraversalNode<BV, Shape>& node,
                BVHModel<BV>& model1, Transform3<typename BV::S>& tf1,
                const Shape& model2, const Transform3<typename BV::S>& tf2,
                const GJKSolver<typename BV::S>* nsolver,
                const CollisionRequest<typename BV::S>& request,
                CollisionResult<typename BV::S>& result, bool use_refit = false,
                bool refit_bottomup = false);

template <typename BV, typename Shape>
void meshShapeCollisionOrientedNodeLeafTesting(
    int b1, int b2, const BVHModel<BV>* model1, const Shape& model2,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2,
    const GJKSolver<typename BV::S>* nsolver, bool enable_statistics,
    int& num_leaf_tests, const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result);

/// @brief Traversal node for mesh and shape, when mesh BVH is one of the
/// oriented node (OBB, RSS, OBBRSS, kIOS)
template <typename Shape>
class MeshShapeCollisionTraversalNodeOBB
    : public MeshShapeCollisionTraversalNode<OBB<typename Shape::S>, Shape> {
 public:
  MeshShapeCollisionTraversalNodeOBB();
  ~MeshShapeCollisionTraversalNodeOBB() override = default;
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;
};

/// @brief Initialize the traversal node for collision between one mesh and one
/// shape, specialized for OBB type
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeOBB<Shape>& node,
                const BVHModel<OBB<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result);

template <typename Shape>
class MeshShapeCollisionTraversalNodeRSS
    : public MeshShapeCollisionTraversalNode<RSS<typename Shape::S>, Shape> {
 public:
  MeshShapeCollisionTraversalNodeRSS();
  ~MeshShapeCollisionTraversalNodeRSS() override = default;
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;
};

/// @brief Initialize the traversal node for collision between one mesh and one
/// shape, specialized for RSS type
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeRSS<Shape>& node,
                const BVHModel<RSS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result);

template <typename Shape>
class MeshShapeCollisionTraversalNodekIOS
    : public MeshShapeCollisionTraversalNode<kIOS<typename Shape::S>, Shape> {
 public:
  MeshShapeCollisionTraversalNodekIOS();
  ~MeshShapeCollisionTraversalNodekIOS() override = default;
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;
};

/// @brief Initialize the traversal node for collision between one mesh and one
///  shape, specialized for kIOS type
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodekIOS<Shape>& node,
                const BVHModel<kIOS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result);

template <typename Shape>
class MeshShapeCollisionTraversalNodeOBBRSS
    : public MeshShapeCollisionTraversalNode<OBBRSS<typename Shape::S>, Shape> {
 public:
  MeshShapeCollisionTraversalNodeOBBRSS();
  ~MeshShapeCollisionTraversalNodeOBBRSS() override = default;
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;
};

/// @brief Initialize the traversal node for collision between one mesh and one
/// shape, specialized for OBBRSS type
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeOBBRSS<Shape>& node,
                const BVHModel<OBBRSS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result);

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/traversal/collision/mesh_shape_collision_traversal_node-inl.h"

#endif
