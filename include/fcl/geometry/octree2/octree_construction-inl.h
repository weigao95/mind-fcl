//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

namespace fcl {
namespace octree2 {

template <typename S>
void Octree<S>::insertVoxelIntoTree(const OctreeVoxel& voxel) {
  std::uint32_t current_node_idx = 0;
  std::uint8_t current_depth = 0;
  while (true) {
    // We are in inner nodes here
    auto& node = inner_nodes_[current_node_idx];
    const bool node_fully_occupied =
        inner_nodes_fully_occupied_[current_node_idx];
    if (node_fully_occupied) return;

    // We would update the child
    assert(!node_fully_occupied);
    const auto child_index = computeChildIndex(voxel, current_depth);
    const bool is_child_leaf = isChildLayerLeafNode(current_depth);
    std::uint32_t child_vector_index = node.children[child_index];
    if (child_vector_index == kInvalidNodeIndex) {
      if (is_child_leaf) {
        // Assign the index
        child_vector_index = static_cast<std::uint32_t>(leaf_nodes_.size());
        node.children[child_index] = child_vector_index;

        // Make a new leaf node
        OctreeLeafNode leaf_node;
        assert(leaf_node.child_occupied.is_all_cleared());

        // Insert into it
        const auto child_of_leaf_node =
            computeChildIndex(voxel, current_depth + 1);
        leaf_node.child_occupied.set_i(child_of_leaf_node);

        // Insert and return
        leaf_nodes_.emplace_back(std::move(leaf_node));
        return;
      } else {
        child_vector_index = static_cast<std::uint32_t>(inner_nodes_.size());
        OctreeInnerNode child_node;
        std::fill(child_node.children.begin(), child_node.children.end(),
                  kInvalidNodeIndex);
        inner_nodes_.emplace_back(std::move(child_node));
        inner_nodes_fully_occupied_.push_back(false);

        // Do NOT use node.children here
        inner_nodes_[current_node_idx].children[child_index] =
            child_vector_index;
      }
    }

    // Move to next inner
    current_depth += 1;
    current_node_idx = child_vector_index;

    // Break if next is leaf
    if (is_child_leaf) break;
  }

  // For leaf nodes
  assert(inner_nodes_.size() == inner_nodes_fully_occupied_.size());
  assert(current_node_idx < leaf_nodes_.size());
  assert(current_depth == meta_info_.leaf_node_depth);
  auto& leaf_node = leaf_nodes_[current_node_idx];
  const auto child_index = computeChildIndex(voxel, meta_info_.leaf_node_depth);
  leaf_node.child_occupied.set_i(child_index);
}

namespace internal {

template <typename S>
struct OctreeInnerNodeAuxiliaryInfoUpdater {
  // Input
  const Octree<S>& tree;
  const std::vector<bool>* inner_nodes_pruned{nullptr};
  const std::vector<OctreeLeafNode>& leaf_nodes;

  // Output
  std::vector<bool>& inner_nodes_full;

  // Constructors
  OctreeInnerNodeAuxiliaryInfoUpdater(
      const Octree<S>& tree_in, const std::vector<bool>* inner_nodes_pruned_in,
      const std::vector<OctreeLeafNode>& leaf_nodes_in,
      std::vector<bool>& inner_nodes_fully_occupied_in)
      : tree(tree_in),
        inner_nodes_pruned(inner_nodes_pruned_in),
        leaf_nodes(leaf_nodes_in),
        inner_nodes_full(inner_nodes_fully_occupied_in) {}

  // Checking the input
  bool sanityCheck() const;

  // Actual func; return whether the node is full
  bool updateRecursive(std::uint32_t inner_node_index, std::uint8_t node_depth);
};

template <typename S>
bool OctreeInnerNodeAuxiliaryInfoUpdater<S>::sanityCheck() const {
  if (tree.leaf_nodes().size() != leaf_nodes.size()) return false;
  if (tree.inner_nodes().size() != inner_nodes_full.size()) return false;
  if (inner_nodes_pruned &&
      inner_nodes_pruned->size() != tree.inner_nodes().size())
    return false;

  // Done
  return true;
}

template <typename S>
bool OctreeInnerNodeAuxiliaryInfoUpdater<S>::updateRecursive(
    std::uint32_t inner_node_index, std::uint8_t node_depth) {
  // Sanity check
  const auto& inner_nodes = tree.inner_nodes();
  assert(inner_nodes.size() == inner_nodes_full.size());
  assert((inner_nodes_pruned == nullptr) ||
         (inner_nodes_pruned->size() == inner_nodes.size()));

  // If pruned
  const OctreeInnerNode& node = inner_nodes[inner_node_index];
  if (inner_nodes_pruned && inner_nodes_pruned->operator[](inner_node_index)) {
    inner_nodes_full[inner_node_index] = false;
    return false;
  }

  // Check nodes
  if (tree.isChildLayerLeafNode(node_depth)) {
    bool is_all_child_occupied = true;
    for (std::uint8_t child_index = 0; child_index < 8; child_index++) {
      const std::uint32_t child_vector_index = node.children[child_index];
      if (child_vector_index == kInvalidNodeIndex) {
        is_all_child_occupied = false;
        break;
      }

      // Else query the child
      const auto child_node = leaf_nodes[child_vector_index];
      const bool child_fully_occupied = child_node.is_fully_occupied();
      if (!child_fully_occupied) {
        is_all_child_occupied = false;
        break;
      }
    }

    // End visiting this node
    inner_nodes_full[inner_node_index] = is_all_child_occupied;
    return is_all_child_occupied;
  }

  // Next layer is also internal
  assert(!tree.isChildLayerLeafNode(node_depth));
  bool is_all_child_occupied = true;
  for (std::uint8_t child_index = 0; child_index < 8; child_index++) {
    const std::uint32_t child_vector_index = node.children[child_index];
    if (child_vector_index == kInvalidNodeIndex) {
      is_all_child_occupied = false;
      continue;  // As we would update its children
    }

    // Call the node recursively
    const bool this_child_node_full =
        updateRecursive(child_vector_index, node_depth + 1);
    if (!this_child_node_full) {
      is_all_child_occupied = false;
    }
  }

  // Done
  inner_nodes_full[inner_node_index] = is_all_child_occupied;
  return is_all_child_occupied;
}

}  // namespace internal

template <typename S>
void Octree<S>::rebuildTree(const PointGenerationFunc& point_generator,
                            int n_points) {
  clearNodes();
  inner_nodes_.reserve(n_points);
  leaf_nodes_.reserve(n_points);
  OctreeVoxel voxel;
  Vector3<S> point;
  S x, y, z;
  for (auto i = 0; i < n_points; i++) {
    // Into voxel
    point_generator(i, x, y, z);
    point.x() = x;
    point.y() = y;
    point.z() = z;
    const bool in_range = computeVoxelCoordinate(point, voxel);
    if (!in_range) continue;

    // Insert into tree
    leaf_points_AABB_ += point;
    insertVoxelIntoTree(voxel);
  }

  // Update the inner status
  updateInnerNodeAuxiliaryInfo(leaf_nodes_, nullptr,
                               inner_nodes_fully_occupied_);
}

template <typename S>
void Octree<S>::updateInnerNodeAuxiliaryInfo(
    const std::vector<OctreeLeafNode>& new_leaf_nodes,
    const std::vector<bool>* inner_nodes_pruned,
    std::vector<bool>& new_inner_nodes_fully_occupied) const {
  // Note: these two might be the same
  if (new_inner_nodes_fully_occupied.size() !=
      inner_nodes_fully_occupied_.size()) {
    new_inner_nodes_fully_occupied = inner_nodes_fully_occupied_;
  }

  // Make updater
  internal::OctreeInnerNodeAuxiliaryInfoUpdater<S> updater(
      *this, inner_nodes_pruned, new_leaf_nodes,
      new_inner_nodes_fully_occupied);
  assert(updater.sanityCheck());
  updater.updateRecursive(0, 0);
}

template <typename S>
bool Octree<S>::Test_insertPointIntoTree(const Vector3<S>& point) {
  // Into voxel
  OctreeVoxel voxel;
  const bool in_range = computeVoxelCoordinate(point, voxel);
  if (!in_range) return false;

  // Make insert
  leaf_points_AABB_ += point;
  insertVoxelIntoTree(voxel);
  return true;
}

template <typename S>
void Octree<S>::rebuildAccordingToPruneInfo(
    const fcl::octree2::OctreePruneInfo& prune_info) {
  rebuildAccordingToPruneInfo(prune_info.prune_internal_nodes,
                              prune_info.new_leaf_nodes);
}

template <typename S>
void Octree<S>::rebuildAccordingToPruneInfo(
    const std::vector<bool>& inner_node_pruned,
    const std::vector<OctreeLeafNode>& leaf_nodes) {
  std::vector<OctreeInnerNode> new_inner_nodes;
  std::vector<OctreeLeafNode> new_leaf_nodes;
  rebuildAccordingToPruneInfo(inner_node_pruned, leaf_nodes, new_inner_nodes,
                              new_leaf_nodes);
  if (new_inner_nodes.empty()) {
    clearNodes();
    return;
  }

  // Assign myself
  inner_nodes_ = std::move(new_inner_nodes);
  leaf_nodes_ = std::move(new_leaf_nodes);
  inner_nodes_fully_occupied_.resize(inner_nodes_.size());
  std::fill(inner_nodes_fully_occupied_.begin(),
            inner_nodes_fully_occupied_.end(), false);
  updateInnerNodeAuxiliaryInfo(leaf_nodes_, nullptr,
                               inner_nodes_fully_occupied_);
}

template <typename S>
void Octree<S>::rebuildAccordingToPruneInfo(
    const std::vector<bool>& inner_node_pruned,
    const std::vector<OctreeLeafNode>& leaf_nodes,
    std::vector<OctreeInnerNode>& new_inner_nodes,
    std::vector<OctreeLeafNode>& new_leaf_nodes) const {
  // The traverse stack
  struct StackElement {
    std::uint32_t node_vector_index;
    std::uint8_t depth;
    std::uint32_t intended_placement;
  };

  // Target output
  new_inner_nodes.clear();
  new_leaf_nodes.clear();
  assert(inner_node_pruned.size() == inner_nodes_.size());
  if (inner_node_pruned[0]) return;

  // Make the stack
  std::stack<StackElement> task_stack;
  StackElement root;
  root.node_vector_index = 0;
  root.depth = 0;
  root.intended_placement = 0;
  task_stack.push(root);
  new_inner_nodes.reserve(inner_nodes_.size());
  new_leaf_nodes.reserve(leaf_nodes.size());
  new_inner_nodes.resize(1);  // For root

  while (!task_stack.empty()) {
    const StackElement this_task = task_stack.top();
    task_stack.pop();

    // Gather info
    const auto node_vector_index = this_task.node_vector_index;
    const auto intended_placement = this_task.intended_placement;
    const bool is_child_leaf = isChildLayerLeafNode(this_task.depth);
    assert(intended_placement < new_inner_nodes.size());
    assert(node_vector_index < inner_nodes_.size());

    // Only handle inner node in the loop
    assert(!inner_node_pruned[node_vector_index]);
    const OctreeInnerNode& node = inner_nodes_[node_vector_index];
    OctreeInnerNode new_inner_node;
    std::fill(new_inner_node.children.begin(), new_inner_node.children.end(),
              kInvalidNodeIndex);

    // Its child is leaf, but it is still inner
    if (is_child_leaf) {
      for (std::uint8_t i = 0; i < 8; i++) {
        const std::uint32_t child_vector_index = node.children[i];
        if (child_vector_index == kInvalidNodeIndex) {
          continue;
        }

        // This is a valid child
        assert(child_vector_index < leaf_nodes.size());
        const auto new_child_vector_index =
            static_cast<std::uint32_t>(new_leaf_nodes.size());
        new_inner_node.children[i] = new_child_vector_index;
        new_leaf_nodes.push_back(leaf_nodes[child_vector_index]);
      }

      // Done with case where child is a leaf
      new_inner_nodes[intended_placement] = new_inner_node;
      continue;
    }

    // This is a leaf
    assert(!is_child_leaf);
    for (std::uint8_t i = 0; i < 8; i++) {
      const std::uint32_t child_vector_index = node.children[i];
      if (child_vector_index == kInvalidNodeIndex) {
        continue;
      }

      // Check is this node pruned
      assert(child_vector_index < inner_node_pruned.size());
      if (inner_node_pruned[child_vector_index]) continue;

      // This is a valid child
      const auto new_child_vector_index =
          static_cast<std::uint32_t>(new_inner_nodes.size());
      new_inner_node.children[i] = new_child_vector_index;
      new_inner_nodes.push_back({});

      // Into stack
      StackElement child_element;
      child_element.node_vector_index = child_vector_index;
      child_element.depth = this_task.depth + 1;
      child_element.intended_placement = new_child_vector_index;
      task_stack.push(std::move(child_element));
    }

    // Assign elem
    new_inner_nodes[intended_placement] = new_inner_node;
  }
}

}  // namespace octree2
}  // namespace fcl
