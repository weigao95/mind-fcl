//
// Created by mech-mind_gw on 3/22/2024.
//

#pragma once

#include <array>
#include <bitset>
#include <limits>

#include "fcl/geometry/octree2/octree_bitset.h"
#include "fcl/geometry/octree2/octree_typedef.h"

namespace fcl {
namespace octree2 {

/// Represent a 2x2x2 unit within the octree, thus it is actually not "leaf"
/// from the perspective of interface. Please refer to OctreeLayerNodeType for
/// detailed description.
struct OctreeLeafNode {
  /// Only contains info about its child, the OctreeLeafNode take a byte
  Bitset8 child_occupied{Bitset8::all_clear};

  // clang-format off
  inline bool is_fully_occupied() const { return child_occupied.is_all_set(); }
  inline bool is_empty() const { return child_occupied.is_all_cleared(); }
  OctreeLayerNodeType node_type() const { return OctreeLayerNodeType::LeafNode; }
  // clang-format on
};

/// The internal node of an octree, the child of an internal node may be
/// InternalNode or LeafNode, depends on current depth of the node.
/// The children are represented as vector index in a flatten vector maintained
/// in the octree.
struct OctreeInnerNode {
  // Children of this node
  using Children = std::array<OctreeNodeIndex, kNumberChildOfOctant>;
  Children children{kInvalidNodeIndex};
};

/// Collection of meta-info for each layer of Octree. This layers are organized
/// from a top-down (coarse to fine) order, with depth (layer index) starts
/// from 0. The bottom layer corresponds to the resolution specified by the
/// user, but no nodes is defined in this layer (RealLeafNoNode).
template <typename S>
struct OctreeLayerMeta {
  std::uint8_t depth{0};
  std::uint16_t full_shape{0};
  std::uint16_t half_shape{0};
  Vector3<S> resolution_xyz{};
  Vector3<S> inv_resolution_xyz{};

  // For nodes
  OctreeLayerNodeType node_type{OctreeLayerNodeType::RealLeafNoNode};
};

/// An octree might be pruned by removing part of the geometry from the tree.
/// This might happen internally in the tree or in this case, as an external
/// modification of the tree.
/// Typically, this struct is used with the original Octree and jointly
/// maintained in, for instance, OctreeCollisionGeometry.
struct OctreePruneInfo {
  // Update on the inner nodes
  // If an inner node is pruned, all its children are ignored
  // Else the inner not is not FULLY pruned (but might be partially pruned)
  // Then its auxiliary info would be updated
  std::vector<bool> prune_internal_nodes{};
  std::vector<bool> new_inner_nodes_fully_occupied{};

  // Replacement of the original leaf nodes
  std::vector<OctreeLeafNode> new_leaf_nodes{};
};

}  // namespace octree2
}  // namespace fcl
