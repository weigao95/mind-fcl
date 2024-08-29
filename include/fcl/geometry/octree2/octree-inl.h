#pragma once

#include <stack>

namespace fcl {
namespace octree2 {

template <typename S>
Octree<S>::Octree(S bottom_resolution_xyz, std::uint16_t bottom_half_shape)
    : Octree(Vector3<S>(bottom_resolution_xyz, bottom_resolution_xyz,
                        bottom_resolution_xyz),
             bottom_half_shape) {}

template <typename S>
Octree<S>::Octree(const Vector3<S>& bottom_resolution_xyz,
                  std::uint16_t bottom_half_shape) {
  std::uint16_t n_half_shape = bottom_half_shape;
  assert(n_half_shape > 0 && ((n_half_shape & (n_half_shape - 1)) == 0));
  assert(n_half_shape >= 2);  // At least three layers

  // As layers_[0] is the top, we must construct in reversed order
  // To make it reversed, compute the top level resolution and shape
  Vector3<S> resolution_xyz = bottom_resolution_xyz;
  while (true) {
    if (n_half_shape > 1) {
      // Update the size and resolution
      resolution_xyz *= 2;
      n_half_shape = n_half_shape / 2;
    } else {
      break;
    }
  }
  assert(n_half_shape == 1);

  // Insert the root layer
  {
    OctreeLayerMeta<S> top_level_meta;
    top_level_meta.depth = 0;
    top_level_meta.full_shape = 1;
    top_level_meta.half_shape = 0;
    top_level_meta.resolution_xyz = resolution_xyz * 2;
    for (auto i = 0; i < 3; i++)
      top_level_meta.inv_resolution_xyz[i] =
          S(1.0) / top_level_meta.resolution_xyz[i];
    top_level_meta.node_type = OctreeLayerNodeType::InnerNode;
    layer_configuration_.emplace_back(std::move(top_level_meta));
  }

  // Insert other layers
  std::uint8_t depth = 1;
  while (true) {
    // Make a new layer
    OctreeLayerMeta<S> this_layer_meta;
    this_layer_meta.depth = depth;
    this_layer_meta.full_shape = n_half_shape * 2;
    this_layer_meta.half_shape = n_half_shape;
    this_layer_meta.resolution_xyz = resolution_xyz;
    this_layer_meta.node_type = OctreeLayerNodeType::InnerNode;
    for (auto i = 0; i < 3; i++)
      this_layer_meta.inv_resolution_xyz[i] =
          S(1.0) / this_layer_meta.resolution_xyz[i];

    // Push into layer configs
    layer_configuration_.emplace_back(std::move(this_layer_meta));

    // Move upward
    if (n_half_shape != bottom_half_shape) {
      resolution_xyz *= S(0.5);
      n_half_shape *= 2;
      depth += 1;
    } else {
      assert(n_half_shape == bottom_half_shape);
      break;
    }
  }

  // Update on node type
  const std::size_t num_layers = layer_configuration_.size();
  assert(num_layers >= 3);
  layer_configuration_[num_layers - 1].node_type =
      OctreeLayerNodeType::RealLeafNoNode;
  layer_configuration_[num_layers - 2].node_type =
      OctreeLayerNodeType::LeafNode;

  // Write the meta info
  meta_info_.num_layers = static_cast<std::uint8_t>(num_layers);
  meta_info_.leaf_node_depth = meta_info_.num_layers - 2;
  meta_info_.real_leaf_resolution = bottom_resolution_xyz;
  meta_info_.leaf_node_resolution = 2 * bottom_resolution_xyz;
  for (auto i = 0; i < 3; i++)
    meta_info_.inv_real_leaf_resolution[i] = S(1.0) / bottom_resolution_xyz[i];
  meta_info_.root_bv.max_ = bottom_resolution_xyz * bottom_half_shape;
  meta_info_.root_bv.min_ = -meta_info_.root_bv.max_;

  // Add root node
  clearNodes();
}

template <typename S>
void Octree<S>::clearNodes() {
  inner_nodes_.clear();
  inner_nodes_fully_occupied_.clear();
  leaf_nodes_.clear();
  leaf_points_AABB_ = AABB<S>();

  // Create an empty root_node
  OctreeInnerNode root_node;
  std::fill(root_node.children.begin(), root_node.children.end(),
            kInvalidNodeIndex);
  inner_nodes_.emplace_back(std::move(root_node));
  inner_nodes_fully_occupied_.push_back(false);
}

template <typename S>
bool Octree<S>::computeVoxelCoordinate(const Vector3<S>& point,
                                       OctreeVoxel& voxel) const {
  const auto& bottom = layer_configuration_.back();
  const auto& inv_resolution = meta_info_.inv_real_leaf_resolution;
  const int layer_half_shape = bottom.half_shape;
  const int layer_full_shape = bottom.full_shape;
  assert(layer_half_shape >= 1);
  using std::floor;
  const int x = floor(point.x() * inv_resolution.x()) + layer_half_shape;
  const int y = floor(point.y() * inv_resolution.y()) + layer_half_shape;
  const int z = floor(point.z() * inv_resolution.z()) + layer_half_shape;
  const bool in_range =
      (x >= 0 && x < layer_full_shape && y >= 0 && y < layer_full_shape &&
       z >= 0 && z < layer_full_shape);
  if (!in_range) return false;

  // Into voxel
  voxel.xyz[0] = static_cast<std::uint16_t>(x);
  voxel.xyz[1] = static_cast<std::uint16_t>(y);
  voxel.xyz[2] = static_cast<std::uint16_t>(z);
  return true;
}

template <typename S>
OctreeChildIndex Octree<S>::computeChildIndex(const Vector3<S>& point,
                                              std::uint8_t parent_layer) const {
  if (parent_layer + 1u >= layer_configuration_.size())
    return kInvalidChildIndex;
  const auto& layer_meta = layer_configuration_[parent_layer + 1];
  const auto& inv_resolution = layer_meta.inv_resolution_xyz;
  const int layer_half_shape = layer_meta.half_shape;
  assert(layer_half_shape >= 1);
  using std::floor;
  const int x = floor(point.x() * inv_resolution.x()) + layer_half_shape;
  const int y = floor(point.y() * inv_resolution.y()) + layer_half_shape;
  const int z = floor(point.z() * inv_resolution.z()) + layer_half_shape;

  // Compute the index
  OctreeChildIndex pos = 0;
  if (x & 1) pos += 1;
  if (y & 1) pos += 2;
  if (z & 1) pos += 4;
  return pos;
}

template <typename S>
bool Octree<S>::isChildLayerLeafNode(std::uint8_t current_depth) const {
  return current_depth + 3 >= meta_info_.num_layers;
}

template <typename S>
OctreeTraverseStackElement<S> Octree<S>::makeStackElementChild(
    const OctreeTraverseStackElement<S>& parent, const AABB<S>& child_AABB,
    std::uint32_t child_vector_index) const {
  OctreeTraverseStackElement<S> child_frame;
  child_frame.bv = child_AABB;
  child_frame.node_vector_index = child_vector_index;
  child_frame.depth = parent.depth + 1;
  child_frame.is_leaf_node = isChildLayerLeafNode(parent.depth);
  return child_frame;
}

template <typename S>
OctreeChildIndex Octree<S>::computeChildIndex(const OctreeVoxel& voxel,
                                              std::uint8_t parent_layer) const {
  const auto num_layers = n_layers();
  if (parent_layer + 1 >= num_layers) return kInvalidChildIndex;
  assert(parent_layer + 1 < num_layers);
  const auto layer_difference = num_layers - parent_layer - 2;
  uint8_t pos = 0;
  if (voxel.x() & (1 << layer_difference)) pos += 1;
  if (voxel.y() & (1 << layer_difference)) pos += 2;
  if (voxel.z() & (1 << layer_difference)) pos += 4;
  return pos;
}

template <typename S>
bool Octree<S>::isPointOccupied(const Vector3<S>& point) const {
  // Into voxel
  OctreeVoxel voxel;
  const bool in_range = computeVoxelCoordinate(point, voxel);
  if (!in_range) return false;

  // Voxel query
  return isVoxelOccupied(voxel);
}

template <typename S>
bool Octree<S>::isVoxelOccupied(const OctreeVoxel& voxel) const {
  // Into as inner nodes
  std::uint32_t current_node_vector_idx = 0;
  std::uint8_t current_node_depth = 0;
  while (true) {
    const auto& node = inner_nodes_[current_node_vector_idx];
    const bool node_full = inner_nodes_fully_occupied_[current_node_vector_idx];
    if (node_full) return true;

    // Into child
    assert(!node_full);
    const bool child_is_leaf_node = isChildLayerLeafNode(current_node_depth);
    const auto child_index = computeChildIndex(voxel, current_node_depth);
    current_node_vector_idx = node.children[child_index];
    current_node_depth += 1;
    assert(current_node_vector_idx != kInvalidNodeIndex);
    if (child_is_leaf_node) break;
  }

  // Must be leaf nodes
  assert(current_node_vector_idx < leaf_nodes_.size());
  const auto& node = leaf_nodes_[current_node_vector_idx];
  const auto child_index = computeChildIndex(voxel, meta_info_.leaf_node_depth);
  return node.child_occupied.test_i(child_index);
}

}  // namespace octree2
}  // namespace fcl
