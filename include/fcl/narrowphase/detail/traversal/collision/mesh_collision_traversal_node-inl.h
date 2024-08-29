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

#ifndef FCL_TRAVERSAL_MESHCOLLISIONTRAVERSALNODE_INL_H
#define FCL_TRAVERSAL_MESHCOLLISIONTRAVERSALNODE_INL_H

#include "fcl/common/unused.h"
#include "fcl/narrowphase/collision_result.h"
#include "fcl/narrowphase/detail/shape_pair_intersect.h"
#include "fcl/narrowphase/detail/traversal/collision/mesh_collision_traversal_node.h"

namespace fcl {

namespace detail {

//==============================================================================
template <typename BV>
MeshCollisionTraversalNode<BV>::MeshCollisionTraversalNode()
    : BVHCollisionTraversalNode<BV>(), nsolver(nullptr) {}

//==============================================================================
template <typename BV>
void MeshCollisionTraversalNode<BV>::leafTesting(int b1, int b2) const {
  if (this->enable_statistics) this->num_leaf_tests++;

  const BVNode<BV>& node1 = this->model1->getBV(b1);
  const BVNode<BV>& node2 = this->model2->getBV(b2);
  const int primitive_id1 = node1.primitiveId();
  const int primitive_id2 = node2.primitiveId();
  const Simplex<S>& simplex1 = this->model1->getSimplex(primitive_id1);
  const Simplex<S>& simplex2 = this->model2->getSimplex(primitive_id2);

  // Make the shape solver
  if (this->model1->isOccupied() && this->model2->isOccupied()) {
    ContactMeta<S> contact;
    contact.o1 = this->model1;
    contact.o2 = this->model2;
    contact.b1 = primitive_id1;
    contact.b2 = primitive_id2;
    ShapePairIntersectSolver<S> shape_solver(nsolver);
    shape_solver.SimplexIntersect(simplex1, Transform3<S>::Identity(), simplex2,
                                  Transform3<S>::Identity(),
                                  Matrix3<S>::Identity(), Vector3<S>::Zero(),
                                  this->request, contact, *(this->result));
  }
}

//==============================================================================
template <typename BV>
bool MeshCollisionTraversalNode<BV>::canStop() const {
  return this->request.terminationConditionSatisfied(*(this->result));
}

//==============================================================================
template <typename BV>
bool initialize(MeshCollisionTraversalNode<BV>& node, BVHModel<BV>& model1,
                Transform3<typename BV::S>& tf1, BVHModel<BV>& model2,
                Transform3<typename BV::S>& tf2,
                const GJKSolver<typename BV::S>* nsolver,
                const CollisionRequest<typename BV::S>& request,
                CollisionResult<typename BV::S>& result, bool use_refit,
                bool refit_bottomup) {
  using S = typename BV::S;

  if (model1.getModelType() != BVH_MODEL_MESH ||
      model2.getModelType() != BVH_MODEL_MESH)
    return false;

  if (!tf1.matrix().isIdentity()) {
    std::vector<Vector3<S>> vertices_transformed1(model1.num_vertices());
    for (int i = 0; i < model1.num_vertices(); ++i) {
      const Vector3<S>& p = model1.getVertex(i);
      Vector3<S> new_v = tf1 * p;
      vertices_transformed1[i] = new_v;
    }

    model1.beginReplaceModel();
    model1.replaceSubModel(vertices_transformed1);
    model1.endReplaceModel(use_refit, refit_bottomup);
    tf1.setIdentity();
  }

  if (!tf2.matrix().isIdentity()) {
    std::vector<Vector3<S>> vertices_transformed2(model2.num_vertices());
    for (int i = 0; i < model2.num_vertices(); ++i) {
      const Vector3<S>& p = model2.getVertex(i);
      Vector3<S> new_v = tf2 * p;
      vertices_transformed2[i] = new_v;
    }

    model2.beginReplaceModel();
    model2.replaceSubModel(vertices_transformed2);
    model2.endReplaceModel(use_refit, refit_bottomup);
    tf2.setIdentity();
  }

  node.model1 = &model1;
  node.tf1 = tf1;
  node.model2 = &model2;
  node.tf2 = tf2;

  node.nsolver = nsolver;
  node.request = request;
  node.result = &result;
  node.cost_density = model1.cost_density * model2.cost_density;
  return true;
}

//==============================================================================
template <typename S>
MeshCollisionTraversalNodeOBB<S>::MeshCollisionTraversalNodeOBB()
    : MeshCollisionTraversalNode<OBB<S>>(), R(Matrix3<S>::Identity()) {
  // Do nothing
}

//==============================================================================
template <typename S>
bool MeshCollisionTraversalNodeOBB<S>::BVTesting(int b1, int b2) const {
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(R, T, this->model1->getBV(b1).bv, this->model2->getBV(b2).bv);
}

//==============================================================================
template <typename S>
void MeshCollisionTraversalNodeOBB<S>::leafTesting(int b1, int b2) const {
  detail::meshCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, this->model2, R, T, this->tf1, this->tf2,
      this->enable_statistics, this->num_leaf_tests, this->nsolver,
      this->request, *this->result);
}

//==============================================================================
template <typename S>
MeshCollisionTraversalNodeRSS<S>::MeshCollisionTraversalNodeRSS()
    : MeshCollisionTraversalNode<RSS<S>>(), R(Matrix3<S>::Identity()) {
  // Do nothing
}

//==============================================================================
template <typename S>
bool MeshCollisionTraversalNodeRSS<S>::BVTesting(int b1, int b2) const {
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(R, T, this->model1->getBV(b1).bv, this->model2->getBV(b2).bv);
}

//==============================================================================
template <typename S>
void MeshCollisionTraversalNodeRSS<S>::leafTesting(int b1, int b2) const {
  detail::meshCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, this->model2, R, T, this->tf1, this->tf2,
      this->enable_statistics, this->num_leaf_tests, this->nsolver,
      this->request, *this->result);
}

//==============================================================================
template <typename S>
MeshCollisionTraversalNodekIOS<S>::MeshCollisionTraversalNodekIOS()
    : MeshCollisionTraversalNode<kIOS<S>>(), R(Matrix3<S>::Identity()) {
  // Do nothing
}

//==============================================================================
template <typename S>
bool MeshCollisionTraversalNodekIOS<S>::BVTesting(int b1, int b2) const {
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(R, T, this->model1->getBV(b1).bv, this->model2->getBV(b2).bv);
}

//==============================================================================
template <typename S>
void MeshCollisionTraversalNodekIOS<S>::leafTesting(int b1, int b2) const {
  detail::meshCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, this->model2, R, T, this->tf1, this->tf2,
      this->enable_statistics, this->num_leaf_tests, this->nsolver,
      this->request, *this->result);
}

//==============================================================================
template <typename S>
MeshCollisionTraversalNodeOBBRSS<S>::MeshCollisionTraversalNodeOBBRSS()
    : MeshCollisionTraversalNode<OBBRSS<S>>(), R(Matrix3<S>::Identity()) {
  // Do nothing
}

//==============================================================================
template <typename S>
bool MeshCollisionTraversalNodeOBBRSS<S>::BVTesting(int b1, int b2) const {
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(R, T, this->model1->getBV(b1).bv, this->model2->getBV(b2).bv);
}

//==============================================================================
template <typename S>
void MeshCollisionTraversalNodeOBBRSS<S>::leafTesting(int b1, int b2) const {
  detail::meshCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, this->model2, R, T, this->tf1, this->tf2,
      this->enable_statistics, this->num_leaf_tests, this->nsolver,
      this->request, *this->result);
}

template <typename BV>
void meshCollisionOrientedNodeLeafTesting(
    int b1, int b2, const BVHModel<BV>* model1, const BVHModel<BV>* model2,
    const Matrix3<typename BV::S>& R, const Vector3<typename BV::S>& T,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2, bool enable_statistics,
    int& num_leaf_tests, const GJKSolver<typename BV::S>* nsolver,
    const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result) {
  using S = typename BV::S;
  (void)(tf2);
  if (enable_statistics) num_leaf_tests++;

  const BVNode<BV>& node1 = model1->getBV(b1);
  const BVNode<BV>& node2 = model2->getBV(b2);
  const int primitive_id1 = node1.primitiveId();
  const int primitive_id2 = node2.primitiveId();
  const Simplex<S>& simplex1 = model1->getSimplex(primitive_id1);
  const Simplex<S>& simplex2 = model2->getSimplex(primitive_id2);

  // Make the shape solver
  if (model1->isOccupied() && model2->isOccupied()) {
    ContactMeta<S> contact;
    contact.o1 = model1;
    contact.o2 = model2;
    contact.b1 = primitive_id1;
    contact.b2 = primitive_id2;

    ShapePairIntersectSolver<S> shape_solver(nsolver);
    shape_solver.SimplexIntersect(simplex1, tf1, simplex2, tf2, R, T, request,
                                  contact, result);
  }
}

//==============================================================================
template <typename BV, typename OrientedNode>
bool setupMeshCollisionOrientedNode(
    OrientedNode& node, const BVHModel<BV>& model1,
    const Transform3<typename BV::S>& tf1, const BVHModel<BV>& model2,
    const Transform3<typename BV::S>& tf2,
    const GJKSolver<typename BV::S>* nsolver,
    const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result) {
  if (model1.getModelType() != BVH_MODEL_MESH ||
      model2.getModelType() != BVH_MODEL_MESH)
    return false;

  node.model1 = &model1;
  node.tf1 = tf1;
  node.model2 = &model2;
  node.tf2 = tf2;

  node.request = request;
  node.result = &result;

  node.nsolver = nsolver;
  node.cost_density = model1.cost_density * model2.cost_density;
  relativeTransform(tf1.linear(), tf1.translation(), tf2.linear(),
                    tf2.translation(), node.R, node.T);
  return true;
}

//==============================================================================
template <typename S>
bool initialize(MeshCollisionTraversalNodeOBB<S>& node,
                const BVHModel<OBB<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<OBB<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result) {
  return detail::setupMeshCollisionOrientedNode(node, model1, tf1, model2, tf2,
                                                nsolver, request, result);
}

//==============================================================================
template <typename S>
bool initialize(MeshCollisionTraversalNodeRSS<S>& node,
                const BVHModel<RSS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<RSS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result) {
  return detail::setupMeshCollisionOrientedNode(node, model1, tf1, model2, tf2,
                                                nsolver, request, result);
}

//==============================================================================
template <typename S>
bool initialize(MeshCollisionTraversalNodekIOS<S>& node,
                const BVHModel<kIOS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<kIOS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result) {
  return detail::setupMeshCollisionOrientedNode(node, model1, tf1, model2, tf2,
                                                nsolver, request, result);
}

//==============================================================================
template <typename S>
bool initialize(MeshCollisionTraversalNodeOBBRSS<S>& node,
                const BVHModel<OBBRSS<S>>& model1, const Transform3<S>& tf1,
                const BVHModel<OBBRSS<S>>& model2, const Transform3<S>& tf2,
                const GJKSolver<S>* nsolver, const CollisionRequest<S>& request,
                CollisionResult<S>& result) {
  return detail::setupMeshCollisionOrientedNode(node, model1, tf1, model2, tf2,
                                                nsolver, request, result);
}

}  // namespace detail
}  // namespace fcl

#endif
