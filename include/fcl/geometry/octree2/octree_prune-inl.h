//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

namespace fcl {
namespace octree2 {

template <typename S>
void pruneOctreeByOBB(const Octree<S>& tree, const OBB<S>& pruned_obb,
                      OctreePruneInfo& existing_prune_info) {
  // Collect the info and init
  const auto& inner_nodes = tree.inner_nodes();
  const auto& original_leaf_nodes = tree.leaf_nodes();
  if (existing_prune_info.prune_internal_nodes.size() == inner_nodes.size()) {
    assert(original_leaf_nodes.size() ==
           existing_prune_info.new_leaf_nodes.size());
    assert(existing_prune_info.new_inner_nodes_fully_occupied.size() ==
           inner_nodes.size());
  } else {
    // Construct a new one
    existing_prune_info.prune_internal_nodes.resize(inner_nodes.size());
    std::fill(existing_prune_info.prune_internal_nodes.begin(),
              existing_prune_info.prune_internal_nodes.end(), false);
    existing_prune_info.new_inner_nodes_fully_occupied =
        tree.inner_nodes_fully_occupied();
    existing_prune_info.new_leaf_nodes = original_leaf_nodes;
  }

  // Input for stacked processing
  auto& leaf_nodes = existing_prune_info.new_leaf_nodes;
  auto& prune_internal_nodes = existing_prune_info.prune_internal_nodes;

  // Make the stack
  std::stack<OctreeTraverseStackElement<S>> task_stack;
  task_stack.push(OctreeTraverseStackElement<S>::MakeRoot(tree.root_bv()));

  // Process loop
  AABB<S> local_aabb;
  OBB<S> obb_for_octree_node;
  obb_for_octree_node.axis.setIdentity();
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task = task_stack.top();
    task_stack.pop();

    // Already pruned
    if ((!this_task.is_leaf_node) &&
        prune_internal_nodes[this_task.node_vector_index]) {
      continue;
    }

    // Make bv for this task and check if this is not related
    obb_for_octree_node.To = this_task.bv.center();
    obb_for_octree_node.extent =
        S(0.5) * (this_task.bv.max_ - this_task.bv.min_);
    if (!pruned_obb.overlap(obb_for_octree_node)) {
      // No overlap, would not prune anything below
      continue;
    }

    // Leaf condition
    if (this_task.is_leaf_node) {
      OctreeLeafNode& leaf_node = leaf_nodes[this_task.node_vector_index];
      for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
        if (!leaf_node.child_occupied.test_i(child_i)) continue;
        computeChildAABB(this_task.bv, child_i, local_aabb);
        const bool to_prune = pruned_obb.contain(local_aabb.center());
        if (to_prune) leaf_node.child_occupied.clear_i(child_i);
      }

      // End for leaf nodes
      continue;
    }

    // Complete pruned
    assert(!this_task.is_leaf_node);
    if (is_contained(pruned_obb, this_task.bv)) {
      prune_internal_nodes[this_task.node_vector_index] = true;
      continue;
    }

    // Not fully pruned
    const OctreeInnerNode& node = inner_nodes[this_task.node_vector_index];
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      // Query the child
      const auto child_vector_index = node.children[child_i];
      if (child_vector_index == kInvalidNodeIndex) continue;

      // Push into stack
      computeChildAABB(this_task.bv, child_i, local_aabb);
      auto child_frame =
          tree.makeStackElementChild(this_task, local_aabb, child_vector_index);
      task_stack.push(std::move(child_frame));
    }
  }

  // Rebuild the meta
  auto& inner_nodes_full = existing_prune_info.new_inner_nodes_fully_occupied;
  tree.updateInnerNodeAuxiliaryInfo(leaf_nodes, &prune_internal_nodes,
                                    inner_nodes_full);
}

}  // namespace octree2
}  // namespace fcl
