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

#ifndef FCL_BVH_MODEL_INL_H
#define FCL_BVH_MODEL_INL_H

#include <algorithm>
#include <new>
#include <queue>
#include <stack>

#include "fcl/geometry/bvh/BVH_model.h"

namespace fcl {

//==============================================================================
template <typename S>
OBJECT_TYPE BVHModelBase<S>::getObjectType() const {
  return OT_BVH;
}

//==============================================================================
template <typename BV>
BVHModelType BVHModel<BV>::getModelType() const {
  if (num_vertices() <= 0) {
    return BVH_MODEL_UNKNOWN;
  }

  // Must be valid
  return num_simplex() > 0 ? BVH_MODEL_MESH : BVH_MODEL_UNKNOWN;
}

//==============================================================================
template <typename BV>
BVHModel<BV>::BVHModel()
    : build_state_(BVH_BUILD_STATE_EMPTY),
      bv_splitter_(new detail::BVSplitter<BV>(detail::SPLIT_METHOD_MEAN)),
      bv_fitter_(new detail::BVFitter<BV>()),
      num_vertex_updated_(0) {
  // Do nothing
  assert(vertex_vector_.empty());
  assert(simplex_vector_.empty());
}

//==============================================================================
template <typename BV>
BVHModel<BV>::BVHModel(const BVHModel<BV>& other)
    : BVHModelBase<S>(other),
      vertex_vector_(other.vertex_vector_),
      simplex_vector_(other.simplex_vector_),
      build_state_(other.build_state_),
      bv_splitter_(other.bv_splitter_),
      bv_fitter_(other.bv_fitter_),
      num_vertex_updated_(other.num_vertex_updated_),
      primitive_indices_(other.primitive_indices_),
      bvs_(other.bvs_) {
  // Do nothing
}

//==============================================================================
template <typename BV>
BVHModel<BV>::~BVHModel() {}

//==============================================================================
template <typename BV>
const BVNode<BV>& BVHModel<BV>::getBV(int id) const {
  return bvs_[id];
}

//==============================================================================
template <typename BV>
BVNode<BV>& BVHModel<BV>::getBV(int id) {
  return bvs_[id];
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::getNumBVs() const {
  return static_cast<int>(bvs_.size());
}

//==============================================================================
template <typename BV>
struct GetNodeTypeImpl {
  static NODE_TYPE run() { return BV_UNKNOWN; }
};

//==============================================================================
template <typename BV>
NODE_TYPE BVHModel<BV>::getNodeType() const {
  return GetNodeTypeImpl<BV>::run();
}

//==============================================================================
template <typename BV>
void BVHModel<BV>::resetSplitMethodForBV(detail::SplitMethodType split_method) {
  bv_splitter_ = std::make_shared<detail::BVSplitter<BV>>(split_method);
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::beginModel(int num_simplex_hint, int num_vertices_hint) {
  if (build_state_ != BVH_BUILD_STATE_EMPTY) {
    // Clear the state
    vertex_vector_.clear();
    simplex_vector_.clear();
    bvs_.clear();
    primitive_indices_.clear();
  }

  // Update on hint size
  if (num_simplex_hint < 0) num_simplex_hint = 8;
  if (num_vertices_hint <= 0) num_vertices_hint = 8;

  // Reserve buffer
  vertex_vector_.reserve(num_vertices_hint);
  if (num_simplex_hint > 0) {
    simplex_vector_.reserve(num_simplex_hint);
  }

  if (build_state_ != BVH_BUILD_STATE_EMPTY) {
    std::cerr << "BVH Warning! Call beginModel() on a BVHModel that is not "
                 "empty. This model was cleared and previous "
                 "triangles/vertices_ were lost.\n";
    build_state_ = BVH_BUILD_STATE_EMPTY;
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Done
  build_state_ = BVH_BUILD_STATE_BEGUN;
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::addVertex(const Vector3<S>& p) {
  if (build_state_ != BVH_BUILD_STATE_BEGUN) {
    std::cerr << "BVH Warning! Call addVertex() in a wrong order. addVertex() "
                 "was ignored. Must do a beginModel() to clear the model for "
                 "addition of new vertices_.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Add point
  vertex_vector_.push_back(p);
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::addTriangle(const Vector3<S>& p1, const Vector3<S>& p2,
                              const Vector3<S>& p3) {
  if (build_state_ == BVH_BUILD_STATE_PROCESSED) {
    std::cerr << "BVH Warning! Call addTriangle() in a wrong order. "
                 "addTriangle() was ignored. Must do a beginModel() to clear "
                 "the model for addition of new triangles.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Push the vertex
  const int offset = static_cast<int>(vertex_vector_.size());
  vertex_vector_.push_back(p1);
  vertex_vector_.push_back(p2);
  vertex_vector_.push_back(p3);

  // Push the new triangle
  MeshSimplex new_triangle(offset, offset + 1, offset + 2);
  simplex_vector_.emplace_back(new_triangle);
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::addTetrahedron(const Vector3<S>& p1, const Vector3<S>& p2,
                                 const Vector3<S>& p3, const Vector3<S>& p4) {
  if (build_state_ == BVH_BUILD_STATE_PROCESSED) {
    std::cerr << "BVH Warning! Call addTriangle() in a wrong order. "
                 "addTriangle() was ignored. Must do a beginModel() to clear "
                 "the model for addition of new triangles.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Push the vertex
  const int offset = static_cast<int>(vertex_vector_.size());
  vertex_vector_.push_back(p1);
  vertex_vector_.push_back(p2);
  vertex_vector_.push_back(p3);
  vertex_vector_.push_back(p4);

  // Push the new triangle
  MeshSimplex new_tetrahedron(offset, offset + 1, offset + 2, offset + 3);
  simplex_vector_.emplace_back(new_tetrahedron);
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::addSubModel(const std::vector<Vector3<S>>& point_vec) {
  if (build_state_ == BVH_BUILD_STATE_PROCESSED) {
    std::cerr << "BVH Warning! Call addSubModel() in a wrong order. "
                 "addSubModel() was ignored. Must do a beginModel() to clear "
                 "the model for addition of new vertices_.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Add points
  for (std::size_t i = 0; i < point_vec.size(); ++i) {
    vertex_vector_.push_back(point_vec[i]);
  }

  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::addSubModel(const std::vector<Vector3<S>>& point_vec,
                              const std::vector<MeshSimplex>& simplex_vec) {
  if (build_state_ == BVH_BUILD_STATE_PROCESSED) {
    std::cerr << "BVH Warning! Call addSubModel() in a wrong order. "
                 "addSubModel() was ignored. Must do a beginModel() to clear "
                 "the model for addition of new vertices_.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // Add points
  const int offset = static_cast<int>(vertex_vector_.size());
  for (std::size_t i = 0; i < point_vec.size(); ++i) {
    vertex_vector_.push_back(point_vec[i]);
  }

  // Add simplex
  for (std::size_t i = 0; i < simplex_vec.size(); ++i) {
    const MeshSimplex& simplex_i = simplex_vec[i];
    MeshSimplex new_simplex =
        simplex_i.is_triangle()
            ? MeshSimplex(simplex_i[0] + offset, simplex_i[1] + offset,
                          simplex_i[2] + offset)
            : MeshSimplex(simplex_i[0] + offset, simplex_i[1] + offset,
                          simplex_i[2] + offset, simplex_i[3] + offset);
    simplex_vector_.emplace_back(new_simplex);
  }

  // Done
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::endModel() {
  if (build_state_ != BVH_BUILD_STATE_BEGUN) {
    std::cerr << "BVH Warning! Call endModel() in wrong order. endModel() was "
                 "ignored.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  if (simplex_vector_.empty() && vertex_vector_.empty()) {
    std::cerr << "BVH Error! endModel() called on model with no triangles and "
                 "vertices_.\n";
    return BVH_ERR_BUILD_EMPTY_MODEL;
  }

  // Build tree from stretch
  buildTree();

  // finish constructing
  build_state_ = BVH_BUILD_STATE_PROCESSED;
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::beginReplaceModel() {
  if (build_state_ != BVH_BUILD_STATE_PROCESSED) {
    std::cerr << "BVH Error! Call beginReplaceModel() on a BVHModel that has "
                 "no previous frame.\n";
    return BVH_ERR_BUILD_EMPTY_PREVIOUS_FRAME;
  }

  num_vertex_updated_ = 0;
  build_state_ = BVH_BUILD_STATE_REPLACE_BEGUN;
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::replaceVertex(const Vector3<S>& p) {
  if (build_state_ != BVH_BUILD_STATE_REPLACE_BEGUN) {
    std::cerr << "BVH Warning! Call replaceVertex() in a wrong order. "
                 "replaceVertex() was ignored. Must do a beginReplaceModel() "
                 "for initialization.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  vertex_vector_[num_vertex_updated_] = p;
  num_vertex_updated_++;
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::replaceSubModel(const std::vector<Vector3<S>>& point_vec) {
  if (build_state_ != BVH_BUILD_STATE_REPLACE_BEGUN) {
    std::cerr << "BVH Warning! Call replaceSubModel() in a wrong order. "
                 "replaceSubModel() was ignored. Must do a beginReplaceModel() "
                 "for initialization.\n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  for (unsigned int i = 0; i < point_vec.size(); ++i) {
    vertex_vector_[num_vertex_updated_] = point_vec[i];
    num_vertex_updated_++;
  }
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::endReplaceModel(bool refit, bool bottomup) {
  if (build_state_ != BVH_BUILD_STATE_REPLACE_BEGUN) {
    std::cerr << "BVH Warning! Call endReplaceModel() in a wrong order. "
                 "endReplaceModel() was ignored. \n";
    return BVH_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  if (num_vertex_updated_ != static_cast<int>(vertex_vector_.size())) {
    std::cerr << "BVH Error! The replaced model should have the same number of "
                 "vertices_ as the old model.\n";
    return BVH_ERR_INCORRECT_DATA;
  }

  if (refit) {
    // refit, do not change BVH structure
    refitTree(bottomup);
  } else {
    // reconstruct bvh tree based on current frame data
    buildTree();
  }

  build_state_ = BVH_BUILD_STATE_PROCESSED;
  return BVH_OK;
}

//==============================================================================
template <typename BV>
Simplex<typename BV::S> BVHModel<BV>::getSimplex(int primitive_id) const {
  const auto& indices_simplex = simplex_vector_[primitive_id];
  if (indices_simplex.is_triangle()) {
    const auto& p0 = vertex_vector_[indices_simplex[0]];
    const auto& p1 = vertex_vector_[indices_simplex[1]];
    const auto& p2 = vertex_vector_[indices_simplex[2]];
    return Simplex<typename BV::S>(p0, p1, p2);
  } else {
    const auto& p0 = vertex_vector_[indices_simplex[0]];
    const auto& p1 = vertex_vector_[indices_simplex[1]];
    const auto& p2 = vertex_vector_[indices_simplex[2]];
    const auto& p3 = vertex_vector_[indices_simplex[3]];
    return Simplex<typename BV::S>(p0, p1, p2, p3);
  }
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::buildTree() {
  return buildTreeIterativeImpl();
}

template <typename BV>
int BVHModel<BV>::buildTreeIterativeImpl() {
  // set BVFitter ans SplitRule
  bv_fitter_->set(vertex_vector_.data(), simplex_vector_.data(),
                  getModelType());
  bv_splitter_->set(vertex_vector_.data(), simplex_vector_.data(),
                    getModelType());

  // Reserve the buffer for bvh
  const std::size_t num_bvs_to_be_allocated =
      simplex_vector_.empty() ? (2 * vertex_vector_.size() - 1)
                              : (2 * simplex_vector_.size() - 1);
  bvs_.reserve(num_bvs_to_be_allocated);
  bvs_.clear();

  // Construct buffer primitive index
  int num_primitives = 0;
  switch (getModelType()) {
    case BVH_MODEL_MESH:
      num_primitives = simplex_vector_.size();
      break;
    case BVH_MODEL_POINTCLOUD:
      num_primitives = vertex_vector_.size();
      break;
    default:
      std::cerr << "BVH Error: Model type not supported!\n";
      return BVH_ERR_UNSUPPORTED_FUNCTION;
  }

  // Set for primitive
  primitive_indices_.resize(num_primitives);
  for (int i = 0; i < num_primitives; ++i) {
    primitive_indices_[i] = i;
  }

  // Build the tree
  iterativeBuildTree(num_primitives);

  // Clear bv fitter/splitter
  bv_fitter_->clear();
  bv_splitter_->clear();
  num_vertex_updated_ = 0;

  // Build state is updated externally
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::iterativeBuildTree(int num_primitives_in) {
  struct QueueElement {
    int bv_node_id;
    int first_primitive;
    int num_primitives;
  };
  std::stack<QueueElement> tasks;

  // The root element
  bvs_.clear();
  {
    QueueElement root;
    bvs_.emplace_back(BVNode<BV>{});
    root.bv_node_id = 0;
    root.first_primitive = 0;
    root.num_primitives = num_primitives_in;
    tasks.push(std::move(root));
  }

  // The loop
  const BVHModelType type = getModelType();
  while (!tasks.empty()) {
    const auto task = tasks.top();
    tasks.pop();

    // A new node
    assert(task.first_primitive + task.num_primitives <=
           static_cast<int>(primitive_indices_.size()));
    assert(task.bv_node_id < static_cast<int>(bvs_.size()));
    auto& bv_node = bvs_[task.bv_node_id];
    bv_node.first_primitive = task.first_primitive;
    bv_node.num_primitives = task.num_primitives;

    // constructing BV
    unsigned int* cur_primitive_indices =
        primitive_indices_.data() + task.first_primitive;
    bv_node.bv = bv_fitter_->fit(cur_primitive_indices, task.num_primitives);
    bv_splitter_->computeRule(bv_node.bv, cur_primitive_indices,
                              task.num_primitives);

    if (task.num_primitives == 1) {
      // bvs_[task.bv_node_id].first_child = -((*cur_primitive_indices) + 1);
      const unsigned loaded_idx = *cur_primitive_indices;
      const int signed_idx = static_cast<int>(loaded_idx);
      bvs_[task.bv_node_id].first_child = - (signed_idx + 1);
    } else {
      // First update on spliting
      int c1 = 0;
      for (int i = 0; i < task.num_primitives; ++i) {
        Vector3<S> p;
        if (type == BVH_MODEL_POINTCLOUD) {
          p = vertex_vector_[cur_primitive_indices[i]];
        } else if (type == BVH_MODEL_MESH) {
          const MeshSimplex& simplex =
              simplex_vector_[cur_primitive_indices[i]];
          const bool is_triangle = simplex.is_triangle();
          const Vector3<S>& p1 = vertex_vector_[simplex[0]];
          const Vector3<S>& p2 = vertex_vector_[simplex[1]];
          const Vector3<S>& p3 = vertex_vector_[simplex[2]];

          // Compute the center point
          if (is_triangle) {
            p.noalias() = (p1 + p2 + p3) / 3.0;
          } else {
            const Vector3<S>& p4 = vertex_vector_[simplex[3]];
            p.noalias() = (p1 + p2 + p3 + p4) / 4.0;
          }
        } else {
          std::cerr << "BVH Error: Model type not supported!\n";
          return BVH_ERR_UNSUPPORTED_FUNCTION;
        }

        // loop invariant: up to (but not including) index c1 in group 1,
        // then up to (but not including) index i in group 2
        //
        //  [1] [1] [1] [1] [2] [2] [2] [x] [x] ... [x]
        //                   c1          i
        //
        if (bv_splitter_->apply(p)) {
          // in the right side, do nothing
        } else {
          std::swap(cur_primitive_indices[i], cur_primitive_indices[c1]);
          c1++;
        }
      }

      // Recursive down
      if ((c1 == 0) || (c1 == task.num_primitives))
        c1 = task.num_primitives / 2;
      int num_first_half = c1;

      // Construct left first, but push right first to match the ordering
      QueueElement left;
      QueueElement right;
      {
        left.bv_node_id = bvs_.size();
        // Be careful with malloc here, Do NOT use mutable ref
        bvs_.emplace_back(BVNode<BV>{});
        bvs_[task.bv_node_id].first_child = left.bv_node_id;
        left.first_primitive = task.first_primitive;
        left.num_primitives = num_first_half;
      }

      {
        right.bv_node_id = bvs_.size();
        bvs_.emplace_back(BVNode<BV>{});
        right.first_primitive = task.first_primitive + num_first_half;
        right.num_primitives = task.num_primitives - num_first_half;
      }

      // Push right first to match the old ordering
      tasks.push(std::move(right));
      tasks.push(std::move(left));
    }
  }

  // Done
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::refitTree(bool bottomup) {
  if (bottomup)
    return refitTreeBottomUp();
  else
    return refitTreeTopDown();
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::refitTreeBottomUp() {
  int res = recursiveRefitTreeBottomUp(0);
  return res;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::recursiveRefitTreeBottomUp(int bv_id) {
  BVNode<BV>* bvnode = bvs_.data() + bv_id;
  if (bvnode->isLeaf()) {
    BVHModelType type = getModelType();
    int primitive_id = -(bvnode->first_child + 1);
    if (type == BVH_MODEL_POINTCLOUD) {
      BV bv;
      fit(vertex_vector_.data() + primitive_id, 1, bv);
      bvnode->bv = bv;
    } else if (type == BVH_MODEL_MESH) {
      // Collect the vertices to fit bv
      const MeshSimplex& simplex = simplex_vector_[primitive_id];
      const int n_point = simplex.is_triangle() ? 3 : 4;
      Vector3<S> v_simplex[4];  // 4 points despite it might be triangle
      for (int i = 0; i < n_point; ++i) {
        v_simplex[i] = vertex_vector_[simplex[i]];
      }

      // Fit the bv and save back
      BV bv;
      fit(v_simplex, n_point, bv);
      bvnode->bv = bv;
    } else {
      std::cerr << "BVH Error: Model type not supported!\n";
      return BVH_ERR_UNSUPPORTED_FUNCTION;
    }
  } else {
    recursiveRefitTreeBottomUp(bvnode->leftChild());
    recursiveRefitTreeBottomUp(bvnode->rightChild());
    bvnode->bv = bvs_[bvnode->leftChild()].bv + bvs_[bvnode->rightChild()].bv;
  }

  // Done
  return BVH_OK;
}

//==============================================================================
template <typename BV>
int BVHModel<BV>::refitTreeTopDown() {
  bv_fitter_->set(vertex_vector_.data(), simplex_vector_.data(),
                  getModelType());
  for (int i = 0; i < getNumBVs(); ++i) {
    BV bv = bv_fitter_->fit(primitive_indices_.data() + bvs_[i].first_primitive,
                            bvs_[i].num_primitives);
    bvs_[i].bv = bv;
  }

  bv_fitter_->clear();
  return BVH_OK;
}

//==============================================================================
template <typename BV>
void BVHModel<BV>::computeLocalAABB() {
  AABB<S> aabb_;
  for (std::size_t i = 0; i < vertex_vector_.size(); ++i) {
    aabb_ += vertex_vector_[i];
  }

  this->aabb_center = aabb_.center();
  this->aabb_radius = 0;
  for (std::size_t i = 0; i < vertex_vector_.size(); ++i) {
    S r = (this->aabb_center - vertex_vector_[i]).squaredNorm();
    if (r > this->aabb_radius) this->aabb_radius = r;
  }

  this->aabb_radius = sqrt(this->aabb_radius);
  this->aabb_local = aabb_;
}

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<AABB<S>> {
  static NODE_TYPE run() { return BV_AABB; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<OBB<S>> {
  static NODE_TYPE run() { return BV_OBB; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<RSS<S>> {
  static NODE_TYPE run() { return BV_RSS; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<kIOS<S>> {
  static NODE_TYPE run() { return BV_kIOS; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<OBBRSS<S>> {
  static NODE_TYPE run() { return BV_OBBRSS; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<KDOP<S, 16>> {
  static NODE_TYPE run() { return BV_KDOP16; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<KDOP<S, 18>> {
  static NODE_TYPE run() { return BV_KDOP18; }
};

//==============================================================================
template <typename S>
struct GetNodeTypeImpl<KDOP<S, 24>> {
  static NODE_TYPE run() { return BV_KDOP24; }
};

}  // namespace fcl

#endif
