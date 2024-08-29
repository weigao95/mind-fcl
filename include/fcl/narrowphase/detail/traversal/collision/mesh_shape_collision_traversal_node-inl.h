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

#ifndef FCL_TRAVERSAL_MESHSHAPECOLLISIONTRAVERSALNODE_INL_H
#define FCL_TRAVERSAL_MESHSHAPECOLLISIONTRAVERSALNODE_INL_H

#include "fcl/common/unused.h"
#include "fcl/narrowphase/detail/traversal/collision/mesh_shape_collision_traversal_node.h"

namespace fcl {

namespace detail {

//==============================================================================
template <typename BV, typename Shape>
MeshShapeCollisionTraversalNode<BV, Shape>::MeshShapeCollisionTraversalNode()
    : BVHShapeCollisionTraversalNode<BV, Shape>(), nsolver(nullptr) {}

//==============================================================================
template <typename BV, typename Shape>
void MeshShapeCollisionTraversalNode<BV, Shape>::leafTesting(int b1,
                                                             int b2) const {
  FCL_UNUSED(b2);
  if (this->enable_statistics) this->num_leaf_tests++;

  const BVNode<BV>& node = this->model1->getBV(b1);
  const int primitive_id = node.primitiveId();
  const Simplex<S> simplex = this->model1->getSimplex(primitive_id);

  // Make the shape solver
  ContactMeta<S> contact;
  contact.o1 = this->model1;
  contact.o2 = this->model2;
  contact.b1 = primitive_id;
  contact.b2 = Contact<S>::NONE;

  // We need to reverse the normal as contact.o1 is NOT the same as object
  // in the solver below (model2)
  contact.reverse_normal = true;

  ShapePairIntersectSolver<S> shape_solver(nsolver);
  shape_solver.template ShapeSimplexIntersect<Shape>(
      *(this->model2), this->tf2, simplex, Transform3<S>::Identity(),
      this->request, contact, *(this->result));
}

//==============================================================================
template <typename BV, typename Shape>
bool MeshShapeCollisionTraversalNode<BV, Shape>::canStop() const {
  return this->request.terminationConditionSatisfied(*(this->result));
}

//==============================================================================
template <typename BV, typename Shape>
bool initialize(MeshShapeCollisionTraversalNode<BV, Shape>& node,
                BVHModel<BV>& model1, Transform3<typename BV::S>& tf1,
                const Shape& model2, const Transform3<typename BV::S>& tf2,
                const GJKSolver<typename BV::S>* nsolver,
                const CollisionRequest<typename BV::S>& request,
                CollisionResult<typename BV::S>& result, bool use_refit,
                bool refit_bottomup) {
  using S = typename BV::S;

  if (model1.getModelType() != BVH_MODEL_MESH) return false;

  if (!tf1.matrix().isIdentity()) {
    std::vector<Vector3<S>> vertices_transformed(model1.num_vertices());
    for (int i = 0; i < model1.num_vertices(); ++i) {
      const Vector3<S>& p = model1.getVertex(i);
      Vector3<S> new_v = tf1 * p;
      vertices_transformed[i] = new_v;
    }

    model1.beginReplaceModel();
    model1.replaceSubModel(vertices_transformed);
    model1.endReplaceModel(use_refit, refit_bottomup);

    tf1.setIdentity();
  }

  node.model1 = &model1;
  node.tf1 = tf1;
  node.model2 = &model2;
  node.tf2 = tf2;
  node.nsolver = nsolver;

  computeBV(model2, tf2, node.model2_bv);
  node.request = request;
  node.result = &result;
  node.cost_density = model1.cost_density * model2.cost_density;
  return true;
}

//==============================================================================
template <typename BV, typename Shape>
void meshShapeCollisionOrientedNodeLeafTesting(
    int b1, int b2, const BVHModel<BV>* model1, const Shape& model2,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2,
    const GJKSolver<typename BV::S>* nsolver, bool enable_statistics,
    int& num_leaf_tests, const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result) {
  FCL_UNUSED(b2);
  using S = typename BV::S;
  if (enable_statistics) num_leaf_tests++;

  const BVNode<BV>& node = model1->getBV(b1);
  const int primitive_id = node.primitiveId();
  const Simplex<S> simplex = model1->getSimplex(primitive_id);

  // Make the shape solver
  ContactMeta<S> contact;
  contact.o1 = model1;
  contact.o2 = &model2;
  contact.b1 = primitive_id;
  contact.b2 = Contact<S>::NONE;

  // We need to reverse the normal as contact.o1 is NOT the same as object
  // in the solver below (model2)
  contact.reverse_normal = true;

  ShapePairIntersectSolver<S> shape_solver(nsolver);
  shape_solver.template ShapeSimplexIntersect<Shape>(
      model2, tf2, simplex, tf1, request, contact, result);
}

//==============================================================================
template <typename Shape>
MeshShapeCollisionTraversalNodeOBB<Shape>::MeshShapeCollisionTraversalNodeOBB()
    : MeshShapeCollisionTraversalNode<OBB<typename Shape::S>, Shape>() {}

//==============================================================================
template <typename Shape>
bool MeshShapeCollisionTraversalNodeOBB<Shape>::BVTesting(int b1,
                                                          int b2) const {
  FCL_UNUSED(b2);
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(this->tf1.linear(), this->tf1.translation(), this->model2_bv,
                  this->model1->getBV(b1).bv);
}

//==============================================================================
template <typename Shape>
void MeshShapeCollisionTraversalNodeOBB<Shape>::leafTesting(int b1,
                                                            int b2) const {
  detail::meshShapeCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, *(this->model2), this->tf1, this->tf2,
      this->nsolver, this->enable_statistics, this->num_leaf_tests,
      this->request, *(this->result));
}

//==============================================================================
template <typename Shape>
MeshShapeCollisionTraversalNodeRSS<Shape>::MeshShapeCollisionTraversalNodeRSS()
    : MeshShapeCollisionTraversalNode<RSS<typename Shape::S>, Shape>() {}

//==============================================================================
template <typename Shape>
bool MeshShapeCollisionTraversalNodeRSS<Shape>::BVTesting(int b1,
                                                          int b2) const {
  FCL_UNUSED(b2);
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(this->tf1.linear(), this->tf1.translation(), this->model2_bv,
                  this->model1->getBV(b1).bv);
}

//==============================================================================
template <typename Shape>
void MeshShapeCollisionTraversalNodeRSS<Shape>::leafTesting(int b1,
                                                            int b2) const {
  detail::meshShapeCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, *(this->model2), this->tf1, this->tf2,
      this->nsolver, this->enable_statistics, this->num_leaf_tests,
      this->request, *(this->result));
}

//==============================================================================
template <typename Shape>
MeshShapeCollisionTraversalNodekIOS<
    Shape>::MeshShapeCollisionTraversalNodekIOS()
    : MeshShapeCollisionTraversalNode<kIOS<typename Shape::S>, Shape>() {}

//==============================================================================
template <typename Shape>
bool MeshShapeCollisionTraversalNodekIOS<Shape>::BVTesting(int b1,
                                                           int b2) const {
  FCL_UNUSED(b2);
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(this->tf1.linear(), this->tf1.translation(), this->model2_bv,
                  this->model1->getBV(b1).bv);
}

//==============================================================================
template <typename Shape>
void MeshShapeCollisionTraversalNodekIOS<Shape>::leafTesting(int b1,
                                                             int b2) const {
  detail::meshShapeCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, *(this->model2), this->tf1, this->tf2,
      this->nsolver, this->enable_statistics, this->num_leaf_tests,
      this->request, *(this->result));
}

//==============================================================================
template <typename Shape>
MeshShapeCollisionTraversalNodeOBBRSS<
    Shape>::MeshShapeCollisionTraversalNodeOBBRSS()
    : MeshShapeCollisionTraversalNode<OBBRSS<typename Shape::S>, Shape>() {}

//==============================================================================
template <typename Shape>
bool MeshShapeCollisionTraversalNodeOBBRSS<Shape>::BVTesting(int b1,
                                                             int b2) const {
  FCL_UNUSED(b2);
  if (this->enable_statistics) this->num_bv_tests++;
  return !overlap(this->tf1.linear(), this->tf1.translation(), this->model2_bv,
                  this->model1->getBV(b1).bv);
}

//==============================================================================
template <typename Shape>
void MeshShapeCollisionTraversalNodeOBBRSS<Shape>::leafTesting(int b1,
                                                               int b2) const {
  detail::meshShapeCollisionOrientedNodeLeafTesting(
      b1, b2, this->model1, *(this->model2), this->tf1, this->tf2,
      this->nsolver, this->enable_statistics, this->num_leaf_tests,
      this->request, *(this->result));
}

template <typename BV, typename Shape, template <typename> class OrientedNode>
bool setupMeshShapeCollisionOrientedNode(
    OrientedNode<Shape>& node, const BVHModel<BV>& model1,
    const Transform3<typename BV::S>& tf1, const Shape& model2,
    const Transform3<typename BV::S>& tf2,
    const GJKSolver<typename BV::S>* nsolver,
    const CollisionRequest<typename BV::S>& request,
    CollisionResult<typename BV::S>& result) {
  if (model1.getModelType() != BVH_MODEL_MESH) return false;

  node.model1 = &model1;
  node.tf1 = tf1;
  node.model2 = &model2;
  node.tf2 = tf2;
  node.nsolver = nsolver;

  computeBV(model2, tf2, node.model2_bv);
  node.request = request;
  node.result = &result;
  node.cost_density = model1.cost_density * model2.cost_density;
  return true;
}

//==============================================================================
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeOBB<Shape>& node,
                const BVHModel<OBB<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result) {
  return detail::setupMeshShapeCollisionOrientedNode(
      node, model1, tf1, model2, tf2, nsolver, request, result);
}

//==============================================================================
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeRSS<Shape>& node,
                const BVHModel<RSS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result) {
  return detail::setupMeshShapeCollisionOrientedNode(
      node, model1, tf1, model2, tf2, nsolver, request, result);
}

//==============================================================================
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodekIOS<Shape>& node,
                const BVHModel<kIOS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result) {
  return detail::setupMeshShapeCollisionOrientedNode(
      node, model1, tf1, model2, tf2, nsolver, request, result);
}

//==============================================================================
template <typename Shape>
bool initialize(MeshShapeCollisionTraversalNodeOBBRSS<Shape>& node,
                const BVHModel<OBBRSS<typename Shape::S>>& model1,
                const Transform3<typename Shape::S>& tf1, const Shape& model2,
                const Transform3<typename Shape::S>& tf2,
                const GJKSolver<typename Shape::S>* nsolver,
                const CollisionRequest<typename Shape::S>& request,
                CollisionResult<typename Shape::S>& result) {
  return detail::setupMeshShapeCollisionOrientedNode(
      node, model1, tf1, model2, tf2, nsolver, request, result);
}

}  // namespace detail
}  // namespace fcl

#endif
