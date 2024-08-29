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

#ifndef FCL_BVH_MODEL_H
#define FCL_BVH_MODEL_H

#include <memory>
#include <vector>

#include "fcl/geometry/bvh/BVH_internal.h"
#include "fcl/geometry/bvh/BV_node.h"
#include "fcl/geometry/bvh/detail/BV_fitter.h"
#include "fcl/geometry/bvh/detail/BV_splitter.h"
#include "fcl/geometry/collision_geometry.h"
#include "fcl/math/bv/OBB.h"
#include "fcl/math/bv/kDOP.h"

namespace fcl {

template <typename S>
class BVHModelBase : public CollisionGeometry<S> {
 public:
  OBJECT_TYPE getObjectType() const override;
  virtual const MeshSimplex* getSimplex() const = 0;
  virtual const Vector3<S>* getVertices() const = 0;
};

/// @brief A class describing the bounding hierarchy of a mesh model or a point
/// cloud model (which is viewed as a degraded version of mesh)
template <typename BV>
class BVHModel : public BVHModelBase<typename BV::S> {
 public:
  using S = typename BV::S;

  /// @brief Model type described by the instance
  BVHModelType getModelType() const;

  /// @brief Constructing an empty BVH
  BVHModel();

  /// @brief copy from another BVH
  BVHModel(const BVHModel& other);

  /// @brief deconstruction, delete mesh data related.
  ~BVHModel();

  /// @brief We provide getBV() and getNumBVs() because BVH may be compressed
  /// (in future), so we must provide some flexibility here

  /// @brief Access the bv giving the its index
  const BVNode<BV>& getBV(int id) const;

  /// @brief Access the bv giving the its index
  BVNode<BV>& getBV(int id);

  /// @brief Get the number of bv in the BVH
  int getNumBVs() const;

  /// @brief Get the BV type: default is unknown
  NODE_TYPE getNodeType() const override;

  /// @brief Compute the AABB for the BVH, used for broad-phase collision
  void computeLocalAABB() override;

  /// @brief update on splitter, called before begin model
  void resetSplitMethodForBV(detail::SplitMethodType split_method);

  /// @brief Begin a new BVH model
  int beginModel(int num_simplex_hint = 0, int num_vertices_hint = 0);

  /// @brief Add one point in the new BVH model
  int addVertex(const Vector3<S>& p);

  /// @brief Add one triangle/tetrahedron in the new BVH model
  int addTriangle(const Vector3<S>& p1, const Vector3<S>& p2,
                  const Vector3<S>& p3);
  int addTetrahedron(const Vector3<S>& p1, const Vector3<S>& p2,
                     const Vector3<S>& p3, const Vector3<S>& p4);

  /// @brief Add a set of triangles in the new BVH model
  int addSubModel(const std::vector<Vector3<S>>& point_vec,
                  const std::vector<MeshSimplex>& simplex_vec);

  /// @brief Add a set of points in the new BVH model
  int addSubModel(const std::vector<Vector3<S>>& point_vec);

  /// @brief End BVH model construction, will build the bounding volume
  /// hierarchy
  int endModel();

  /// @brief Replace the geometry information of current frame (i.e. should have
  /// the same mesh topology with the previous frame)
  int beginReplaceModel();

  /// @brief Replace one point in the old BVH model
  int replaceVertex(const Vector3<S>& p);

  /// @brief Replace a set of points in the old BVH model
  int replaceSubModel(const std::vector<Vector3<S>>& point_vec);

  /// @brief End BVH model replacement, will also refit or rebuild the bounding
  /// volume hierarchy
  int endReplaceModel(bool refit = true, bool bottomup = true);

  /// Implement the interface
  Simplex<S> getSimplex(int primitive_id) const;
  inline const Vector3<S>& getVertex(int vertex_index) const {
    return vertex_vector_[vertex_index];
  }

  /// @brief Simple const access of internal data
  // clang-format off
  inline int num_vertices() const { return static_cast<int>(vertex_vector_.size()); }
  inline const Vector3<S>* vertices() const { return vertex_vector_.data(); };
  Vector3<S>* mutable_vertices() { return vertex_vector_.data(); };

  inline const MeshSimplex* simplex_indices() const { return simplex_vector_.data(); }
  inline int num_simplex() const { return static_cast<int>(simplex_vector_.size()); }
  inline BVHBuildState build_state() const { return build_state_; }

  // Internal access of primitive id
  const std::vector<unsigned int>& Test_primitiveIndices() const { return primitive_indices_; }
  // clang-format on

  /// @brief access from base class without template on BV
  const MeshSimplex* getSimplex() const override { return simplex_indices(); }
  const Vector3<S>* getVertices() const override { return vertices(); }

 private:
  // raw data of the geometry
  std::vector<Vector3<S>> vertex_vector_;
  std::vector<MeshSimplex> simplex_vector_;

  // Info for building
  BVHBuildState build_state_{BVH_BUILD_STATE_EMPTY};

  // Processor of bv used in construction
  std::shared_ptr<detail::BVSplitterBase<BV>> bv_splitter_;
  std::shared_ptr<detail::BVFitterBase<BV>> bv_fitter_;

  // for ccd vertex update
  int num_vertex_updated_{0};

  // Re-ordering of the input triangles
  std::vector<unsigned int> primitive_indices_;

  /// @brief Bounding volume hierarchy
  std::vector<BVNode<BV>> bvs_;

  /// @brief Build the bounding volume hierarchy
  int buildTree();
  int buildTreeIterativeImpl();
  int iterativeBuildTree(int num_primitives);

  /// @brief Refit the bounding volume hierarchy
  int refitTree(bool bottomup);

  /// @brief Refit the bounding volume hierarchy in a top-down way (slow but
  /// more compact)
  int refitTreeTopDown();

  /// @brief Refit the bounding volume hierarchy in a bottom-up way (fast but
  /// less compact)
  int refitTreeBottomUp();

  /// @brief Recursive kernel for bottomup refitting
  int recursiveRefitTreeBottomUp(int bv_id);
};

}  // namespace fcl

#include "fcl/geometry/bvh/BVH_model-inl.h"

#endif
