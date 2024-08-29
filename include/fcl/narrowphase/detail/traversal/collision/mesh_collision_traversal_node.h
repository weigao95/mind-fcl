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

#ifndef FCL_TRAVERSAL_MESHCOLLISIONTRAVERSALNODE_H
#define FCL_TRAVERSAL_MESHCOLLISIONTRAVERSALNODE_H

#include "fcl/math/bv/OBB.h"
#include "fcl/math/bv/OBBRSS.h"
#include "fcl/math/bv/RSS.h"
#include "fcl/math/bv/kIOS.h"
#include "fcl/narrowphase/contact.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/traversal/collision/bvh_collision_traversal_node.h"
#include "fcl/narrowphase/detail/traversal/collision/intersect.h"

namespace fcl {

namespace detail {

/// @brief Traversal node for collision between two meshes
///        In its general form, mesh collision node can NOT handle any
///        transformation of the shape (for example, AABB does NOT support
///        bv testing with orientation). As a result, we must build a new BVH
///        using the transformed vertices and use the default bv test.
template <typename BV>
class MeshCollisionTraversalNode : public BVHCollisionTraversalNode<BV> {
 public:
  using S = typename BV::S;
  MeshCollisionTraversalNode();
  ~MeshCollisionTraversalNode() override = default;

  /// Implment the interface
  void leafTesting(int b1, int b2) const override;
  bool canStop() const override;

  // Raw data
  S cost_density;
  const GJKSolver<S>* nsolver;
};

/// @brief Initialize traversal node for collision between two meshes, given the
/// current transforms. As mentioned above, for some BV (such as AABB) bv
/// testing with orientation is NOT supported. As a result, we must build a new
/// BVH using the transformed vertices and use the default bv test.
template <typename BV>
bool initialize(MeshCollisionTraversalNode<BV>& node, BVHModel<BV>& model1,
                Transform3<typename BV::S>& tf1, BVHModel<BV>& model2,
                Transform3<typename BV::S>& tf2,
                const GJKSolver<typename BV::S>* nsolver,
                const CollisionRequest<typename BV::S>& request,
                CollisionResult<typename BV::S>& result, bool use_refit = false,
                bool refit_bottomup = false);

/// @brief Traversal node for collision between two meshes if their underlying
/// BVH node is oriented node (OBB, RSS, OBBRSS, kIOS)
template <typename S>
class MeshCollisionTraversalNodeOBB
    : public MeshCollisionTraversalNode<OBB<S>> {
 public:
  MeshCollisionTraversalNodeOBB();
  ~MeshCollisionTraversalNodeOBB() override = default;

  /// Implement the interface
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;

  Matrix3<S> R;
  Vector3<S> T;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @brief Initialize traversal node for collision between two meshes,
/// specialized for OBB type
template <typename S>
bool initialize(MeshCollisionTraversalNodeOBB<S>& node,
                const BVHModel<OBB<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<OBB<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result);

template <typename S>
class MeshCollisionTraversalNodeRSS
    : public MeshCollisionTraversalNode<RSS<S>> {
 public:
  MeshCollisionTraversalNodeRSS();
  ~MeshCollisionTraversalNodeRSS() override = default;

  /// Implement the interface
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;

  Matrix3<S> R;
  Vector3<S> T;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @brief Initialize traversal node for collision between two meshes,
/// specialized for RSS type
template <typename S>
bool initialize(MeshCollisionTraversalNodeRSS<S>& node,
                const BVHModel<RSS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<RSS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result);

template <typename S>
class MeshCollisionTraversalNodekIOS
    : public MeshCollisionTraversalNode<kIOS<S>> {
 public:
  MeshCollisionTraversalNodekIOS();
  ~MeshCollisionTraversalNodekIOS() override = default;

  // Implement the interface
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;

  Matrix3<S> R;
  Vector3<S> T;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @brief Initialize traversal node for collision between two meshes,
/// specialized for kIOS type
template <typename S>
bool initialize(MeshCollisionTraversalNodekIOS<S>& node,
                const BVHModel<kIOS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<kIOS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result);

template <typename S>
class MeshCollisionTraversalNodeOBBRSS
    : public MeshCollisionTraversalNode<OBBRSS<S>> {
 public:
  MeshCollisionTraversalNodeOBBRSS();
  ~MeshCollisionTraversalNodeOBBRSS() override = default;

  // Implement the interface
  bool BVTesting(int b1, int b2) const override;
  void leafTesting(int b1, int b2) const override;

  Matrix3<S> R;
  Vector3<S> T;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @brief Initialize traversal node for collision between two meshes,
/// specialized for OBBRSS type
template <typename S>
bool initialize(MeshCollisionTraversalNodeOBBRSS<S>& node,
                const BVHModel<OBBRSS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<OBBRSS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result);

template <typename BV>
void meshCollisionOrientedNodeLeafTesting(
    int b1, int b2, const BVHModel<BV>* model1, const BVHModel<BV>* model2,
    const Matrix3<typename BV::S>& R, const Vector3<typename BV::S>& T,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2, bool enable_statistics,
    int& num_leaf_tests, const GJKSolver<typename BV::S>* nsolver,
    const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result);

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/traversal/collision/mesh_collision_traversal_node-inl.h"

#endif
