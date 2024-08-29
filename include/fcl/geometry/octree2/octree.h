//
// Created by mech-mind_gw on 3/22/2024.
//

#pragma once

#include <cassert>
#include <vector>

#include "fcl/common/types.h"
#include "fcl/geometry/octree2/octree_node.h"
#include "fcl/geometry/octree2/octree_util.h"

namespace fcl {
namespace octree2 {

/// Octree represents a 3D geometry consists of many boxes. These
/// The represented ranges is:
///    [-half_range_x, half_range_x] x [-h_y, h_y] x [-h_z, h_z]
/// where: half_range_x/y/z = bottom_half_shape * resolution_x/y/z.
/// The bottom_half_shape must be the power of two. On the other hand, it can
/// be very large as the internal data structure is sparse.
///
/// This Octree does NOT intend to support arbitrary update (insert/remove).
/// Instead, this Octree would be built once and pruned later.
/// \tparam S typically float
template <typename S>
class Octree {
 public:
  Octree(const Vector3<S>& bottom_resolution_xyz,
         std::uint16_t bottom_half_shape);
  Octree(S bottom_resolution_xyz, std::uint16_t bottom_half_shape);
  ~Octree() = default;

  // Query of internal info
  // clang-format off
  const std::vector<OctreeInnerNode>& inner_nodes() const { return inner_nodes_; }
  const std::vector<bool>& inner_nodes_fully_occupied() const { return inner_nodes_fully_occupied_; }
  const std::vector<OctreeLeafNode>& leaf_nodes() const { return leaf_nodes_; }
  const std::vector<OctreeLayerMeta<S>>& layer_metas() const { return layer_configuration_; }
  // clang-format on

  // Query of meta info
  std::uint8_t n_layers() const { return meta_info_.num_layers; };
  const AABB<S>& root_bv() const { return meta_info_.root_bv; }
  const AABB<S>& leaf_points_AABB() const { return leaf_points_AABB_; }
  std::size_t n_inner_nodes() const { return inner_nodes_.size(); }
  std::size_t n_leaf_nodes() const { return leaf_nodes_.size(); }
  inline bool isChildLayerLeafNode(std::uint8_t current_layer_index) const;
  OctreeTraverseStackElement<S> makeStackElementChild(
      const OctreeTraverseStackElement<S>& parent, const AABB<S>& child_AABB,
      std::uint32_t child_vector_index) const;
  const Vector3<S>& bottom_resolution_xyz() const {
    return meta_info_.real_leaf_resolution;
  }

  // Query the voxel and point
  // clang-format off
  bool computeVoxelCoordinate(const Vector3<S>& point, OctreeVoxel& voxel) const;
  OctreeChildIndex computeChildIndex(const Vector3<S>& point, std::uint8_t parent_layer) const;
  OctreeChildIndex computeChildIndex(const OctreeVoxel& voxel, std::uint8_t parent_layer) const;
  bool isPointOccupied(const Vector3<S>& point) const;
  bool isVoxelOccupied(const OctreeVoxel& voxel) const;
  // clang-format on

  // Used for building
  void clearNodes();
  using PointGenerationFunc = std::function<void(int index, S& x, S& y, S& z)>;
  void rebuildTree(const PointGenerationFunc& point_generator, int n_points);
  void updateInnerNodeAuxiliaryInfo(
      const std::vector<OctreeLeafNode>& new_leaf_nodes,
      const std::vector<bool>* inner_nodes_pruned,
      std::vector<bool>& new_inner_nodes_fully_occupied) const;
  void rebuildAccordingToPruneInfo(const OctreePruneInfo& prune_info);

 private:
  // Inner and leaf nodes
  std::vector<OctreeInnerNode> inner_nodes_;
  std::vector<bool> inner_nodes_fully_occupied_;
  std::vector<OctreeLeafNode> leaf_nodes_;
  AABB<S> leaf_points_AABB_; // The AABB expanded by inserted points

  // Meta-info
  std::vector<OctreeLayerMeta<S>> layer_configuration_;
  struct {
    std::uint8_t leaf_node_depth{0};
    std::uint8_t num_layers{0};
    Vector3<S> leaf_node_resolution;
    Vector3<S> real_leaf_resolution;
    Vector3<S> inv_real_leaf_resolution;
    AABB<S> root_bv;
  } meta_info_;

  // Internal utility
  void insertVoxelIntoTree(const OctreeVoxel& key);
  void rebuildAccordingToPruneInfo(
      const std::vector<bool>& inner_node_pruned,
      const std::vector<OctreeLeafNode>& leaf_nodes,
      std::vector<OctreeInnerNode>& new_inner_nodes,
      std::vector<OctreeLeafNode>& new_leaf_nodes) const;
  void rebuildAccordingToPruneInfo(
      const std::vector<bool>& inner_node_pruned,
      const std::vector<OctreeLeafNode>& leaf_nodes);

 public:
  // Used for testing
  bool Test_insertPointIntoTree(const Vector3<S>& point);
};

}  // namespace octree2
}  // namespace fcl

#include "fcl/geometry/octree2/octree-inl.h"
#include "fcl/geometry/octree2/octree_construction-inl.h"
