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
/** @author Sean Curtis (2018) Modify API and correct implementation bugs. */

#ifndef FCL_SHAPE_CONVEX_INL_H
#define FCL_SHAPE_CONVEX_INL_H

#include <map>
#include <set>
#include <utility>

#include "fcl/geometry/shape/convex.h"

namespace fcl {

//==============================================================================
template <typename S>
Convex<S>::Convex(
    const std::shared_ptr<const std::vector<Vector3<S>>>& vertices,
    int num_faces, const std::shared_ptr<const std::vector<int>>& faces,
    bool throw_if_invalid)
    : ShapeBase<S>(),
      vertices_(vertices),
      num_faces_(num_faces),
      faces_(faces),
      find_extreme_via_neighbors_{vertices->size() >
                                  kMinVertCountForEdgeWalking} {
  assert(vertices != nullptr);
  assert(faces != nullptr);
  // Compute an interior point. We're computing the mean point and *not* some
  // alternative such as the centroid or bounding box center.
  Vector3<S> sum = Vector3<S>::Zero();
  for (const auto& vertex : *vertices_) {
    sum += vertex;
  }
  interior_point_ = sum * (S)(1.0 / vertices_->size());
  FindVertexNeighbors();
  ValidateMesh(throw_if_invalid);
}

//==============================================================================
template <typename S>
void Convex<S>::computeLocalAABB() {
  this->aabb_local.min_.setConstant(std::numeric_limits<S>::max());
  this->aabb_local.max_.setConstant(-std::numeric_limits<S>::max());

  for (const auto& v : *vertices_) {
    this->aabb_local += v;
  }

  this->aabb_center = this->aabb_local.center();
  this->aabb_radius = (this->aabb_local.min_ - this->aabb_center).norm();
}

//==============================================================================
template <typename S>
NODE_TYPE Convex<S>::getNodeType() const {
  return GEOM_CONVEX;
}

//==============================================================================
template <typename S>
std::vector<Vector3<S>> Convex<S>::getBoundVertices(
    const Transform3<S>& tf) const {
  std::vector<Vector3<S>> result;
  result.reserve(vertices_->size());

  for (const auto& v : *vertices_) {
    result.push_back(tf * v);
  }

  return result;
}

//==============================================================================
template <typename S>
const Vector3<S>& Convex<S>::findExtremeVertex(const Vector3<S>& v_C) const {
  // TODO(SeanCurtis-TRI): Create an override of this that seeds the search with
  //  the last extremal vertex index (assuming some kind of coherency in the
  //  evaluation sequence).
  if (find_extreme_via_neighbors_) {
    return findExtremeVertexViaNeighbours(v_C);
  } else {
    return findExtremeVertexNaive(v_C);
  }
}

//==============================================================================
template <typename S>
const Vector3<S>& Convex<S>::findExtremeVertexNaive(
    const Vector3<S>& v_C) const {
  const std::vector<Vector3<S>>& vertices = *vertices_;
  const auto extreme_index = findExtremeVertexIndexNaive(v_C);
  return vertices[extreme_index];
}

//==============================================================================
template <typename S>
int Convex<S>::findExtremeVertexIndexNaive(const Vector3<S>& v_C) const {
  const std::vector<Vector3<S>>& vertices = *vertices_;
  int extreme_index = 0;
  S extreme_value = v_C.dot(vertices[extreme_index]);

  // Simple linear search.
  for (int i = 1; i < static_cast<int>(vertices.size()); ++i) {
    S value = v_C.dot(vertices[i]);
    if (value > extreme_value) {
      extreme_index = i;
      extreme_value = value;
    }
  }

  // Done
  return extreme_index;
}

//==============================================================================
template <typename S>
void Convex<S>::initExtremeViaNeighborCache() {
  // +x
  {
    const int unit_x_index = findExtremeVertexIndexNaive(Vector3<S>::UnitX());
    init_direction_vertex_cache_[0] =
        std::make_pair(Vector3<S>::UnitX(), unit_x_index);
  }

  // -x
  {
    Vector3<S> negative_unit_x = -Vector3<S>::UnitX();
    const int neg_unit_x_index = findExtremeVertexIndexNaive(negative_unit_x);
    init_direction_vertex_cache_[1] =
        std::make_pair(std::move(negative_unit_x), neg_unit_x_index);
  }

  // +y
  {
    const int unit_y_index = findExtremeVertexIndexNaive(Vector3<S>::UnitY());
    init_direction_vertex_cache_[2] =
        std::make_pair(Vector3<S>::UnitY(), unit_y_index);
  }

  // -y
  {
    Vector3<S> negative_unit_y = -Vector3<S>::UnitY();
    const int neg_unit_y_index = findExtremeVertexIndexNaive(negative_unit_y);
    init_direction_vertex_cache_[3] =
        std::make_pair(std::move(negative_unit_y), neg_unit_y_index);
  }

  // +z
  {
    const int unit_z_index = findExtremeVertexIndexNaive(Vector3<S>::UnitZ());
    init_direction_vertex_cache_[4] =
        std::make_pair(Vector3<S>::UnitZ(), unit_z_index);
  }

  // -z
  {
    Vector3<S> negative_unit_z = -Vector3<S>::UnitZ();
    const int neg_unit_z_index = findExtremeVertexIndexNaive(negative_unit_z);
    init_direction_vertex_cache_[5] =
        std::make_pair(std::move(negative_unit_z), neg_unit_z_index);
  }
}

//==============================================================================
template <typename S>
const Vector3<S>& Convex<S>::findExtremeVertexViaNeighbours(
    const Vector3<S>& v_C) const {
  // Init checking
  assert(find_extreme_via_neighbors_);

  // First search over the init directions
  int init_vertex_index = -1;
  {
    S max_dot = 0;
    for (const auto& init_direction_index : init_direction_vertex_cache_) {
      const Vector3<S>& this_direction = init_direction_index.first;
      const int this_vertex_idx = init_direction_index.second;
      const S dot_value = this_direction.dot(v_C);
      if (init_vertex_index < 0 || (dot_value > max_dot)) {
        init_vertex_index = this_vertex_idx;
        max_dot = dot_value;
      }
    }

    // Must be positive
    assert(max_dot > 0);
  }

  // Init is ready
  assert(init_vertex_index >= 0);
  const std::vector<Vector3<S>>& vertices = *vertices_;
  int extreme_index = init_vertex_index;
  S extreme_value = v_C.dot(vertices[extreme_index]);

  // Cache for searched variable
  int searched_parent_0 = init_vertex_index;
  int searched_parent_1 = init_vertex_index;

  // Search loop
  bool keep_searching = true;
  while (keep_searching) {
    keep_searching = false;
    const int neighbor_start = neighbors_[extreme_index];
    const int neighbor_count = neighbors_[neighbor_start];
    const int old_extreme_index = extreme_index;

    // Start search the neighbors
    for (int n_index = neighbor_start + 1;
         n_index <= neighbor_start + neighbor_count; ++n_index) {
      const int neighbor_index = neighbors_[n_index];

      // Check visited
      if (neighbor_index == searched_parent_0 ||
          neighbor_index == searched_parent_1) {
        continue;
      }

      const S neighbor_value = v_C.dot(vertices[neighbor_index]);
      /// In the old impl, the following comments are provided:
      // N.B. Testing >= (instead of >) protects us from the (very rare) case
      // where the *starting* vertex is co-planar with all of its neighbors
      // *and* the query direction v_C is perpendicular to that plane. With >
      // we wouldn't walk away from the start vertex. With >= we will walk
      // away and continue around (although it won't necessarily be the
      // fastest walk down).
      //
      /// In the new impl, we can use > here due to the init
      if (neighbor_value > extreme_value) {
        // Update cache
        searched_parent_1 = extreme_index;

        // Update flag
        keep_searching = true;
        extreme_index = neighbor_index;
        extreme_value = neighbor_value;
      }
      // else do nothing
    }

    // Update external cache
    searched_parent_0 = old_extreme_index;
  }

  // Done
  return vertices[extreme_index];
}

//==============================================================================
template <typename S>
void Convex<S>::ValidateMesh(bool throw_on_error) {
  ValidateTopology(throw_on_error);
  // TODO(SeanCurtis-TRI) Implement the missing "all-faces-are-planar" test.
  // TODO(SeanCurtis-TRI) Implement the missing "really-is-convex" test.
}

//==============================================================================
template <typename S>
void Convex<S>::ValidateTopology(bool throw_on_error) {
  // Computing the vertex neighbors is a pre-requisite to determining validity.
  assert(neighbors_.size() > vertices_->size());

  std::stringstream ss;
  ss << "Found errors in the Convex mesh:";

  // To simplify the code, we define an edge as a pair of ints (A, B) such that
  // A < B must be true.
  auto make_edge = [](int v0, int v1) {
    if (v0 > v1) std::swap(v0, v1);
    return std::make_pair(v0, v1);
  };

  bool all_connected = true;
  // Build a map from each unique edge to the _number_ of adjacent faces (see
  // the definition of make_edge for the encoding of an edge).
  std::map<std::pair<int, int>, int> per_edge_face_count;
  // First, pre-populate all the edges found in the vertex neighbor calculation.
  for (int v = 0; v < static_cast<int>(vertices_->size()); ++v) {
    const int neighbor_start = neighbors_[v];
    const int neighbor_count = neighbors_[neighbor_start];
    if (neighbor_count == 0) {
      if (all_connected) {
        ss << "\n Not all vertices are connected.";
        all_connected = false;
      }
      ss << "\n  Vertex " << v << " is not included in any faces.";
    }
    for (int n_index = neighbor_start + 1;
         n_index <= neighbor_start + neighbor_count; ++n_index) {
      const int n = neighbors_[n_index];
      per_edge_face_count[make_edge(v, n)] = 0;
    }
  }

  // To count adjacent faces, we must iterate through the faces; we can't infer
  // it from how many times an edge appears in the neighbor list. In the
  // degenerate case where three faces share an edge, that edge would only
  // appear twice. So, we must explicitly examine each face.
  const std::vector<int>& faces = *faces_;
  int face_index = 0;
  for (int f = 0; f < num_faces_; ++f) {
    const int vertex_count = faces[face_index];
    int prev_v = faces[face_index + vertex_count];
    for (int i = face_index + 1; i <= face_index + vertex_count; ++i) {
      const int v = faces[i];
      ++per_edge_face_count[make_edge(v, prev_v)];
      prev_v = v;
    }
    face_index += vertex_count + 1;
  }

  // Now examine the results.
  bool is_watertight = true;
  for (const auto& key_value_pair : per_edge_face_count) {
    const auto& edge = key_value_pair.first;
    const int count = key_value_pair.second;
    if (count != 2) {
      if (is_watertight) {
        is_watertight = false;
        ss << "\n The mesh is not watertight.";
      }
      ss << "\n  Edge between vertices " << edge.first << " and " << edge.second
         << " is shared by " << count << " faces (should be 2).";
    }
  }
  // We can't trust walking the edges on a mesh that isn't watertight.
  const bool has_error = !(is_watertight && all_connected);
  // Note: find_extreme_via_neighbors_ may already be false because the mesh
  // is too small. Don't indirectly re-enable it just because the mesh doesn't
  // have any errors.
  find_extreme_via_neighbors_ = find_extreme_via_neighbors_ && !has_error;
  if (has_error && throw_on_error) {
    throw std::runtime_error(ss.str());
  }

  // Init the cache for findExtreme
  if (find_extreme_via_neighbors_) {
    initExtremeViaNeighborCache();
  }
}

//==============================================================================
template <typename S>
void Convex<S>::FindVertexNeighbors() {
  // We initially build it using sets. Two faces with a common edge will
  // independently want to report that the edge's vertices are neighbors. So,
  // we rely on the set to eliminate the redundant declaration and then dump
  // the results to a more compact representation: the vector.
  const int v_count = static_cast<int>(vertices_->size());
  std::vector<std::set<int>> neighbors(v_count);
  const std::vector<int>& faces = *faces_;
  int face_index = 0;
  for (int f = 0; f < num_faces_; ++f) {
    const int vertex_count = faces[face_index];
    int prev_v = faces[face_index + vertex_count];
    for (int i = face_index + 1; i <= face_index + vertex_count; ++i) {
      const int v = faces[i];
      neighbors[v].insert(prev_v);
      neighbors[prev_v].insert(v);
      prev_v = v;
    }
    face_index += vertex_count + 1;
  }

  // Now build the encoded adjacency graph as documented.
  neighbors_.resize(v_count);
  for (int v = 0; v < v_count; ++v) {
    const std::set<int>& v_neighbors = neighbors[v];
    neighbors_[v] = static_cast<int>(neighbors_.size());
    neighbors_.push_back(static_cast<int>(v_neighbors.size()));
    neighbors_.insert(neighbors_.end(), v_neighbors.begin(), v_neighbors.end());
  }
}

}  // namespace fcl

#endif
