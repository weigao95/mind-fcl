//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

namespace fcl {
namespace octree2 {

template <typename S>
void visitOctree(const Octree<S>& tree,
                 const OctreePruneInfo* prune_octree_info,
                 const VisitOctreeNodeFunc<S>& visitor) {
  // Collection inputs
  const auto& inner_nodes = tree.inner_nodes();
  const auto& inner_nodes_full =
      (prune_octree_info != nullptr &&
       prune_octree_info->new_inner_nodes_fully_occupied.size() ==
           tree.inner_nodes_fully_occupied().size())
          ? prune_octree_info->new_inner_nodes_fully_occupied
          : tree.inner_nodes_fully_occupied();
  const auto& leaf_nodes =
      (prune_octree_info != nullptr &&
       prune_octree_info->new_leaf_nodes.size() == tree.leaf_nodes().size())
          ? prune_octree_info->new_leaf_nodes
          : tree.leaf_nodes();
  const bool try_prune_inner_nodes =
      (prune_octree_info != nullptr) &&
      (prune_octree_info->prune_internal_nodes.size() == inner_nodes.size());
  assert(inner_nodes.size() == inner_nodes_full.size());

  // Make the stack
  std::stack<OctreeTraverseStackElement<S>> task_stack;
  task_stack.push(OctreeTraverseStackElement<S>::MakeRoot(tree.root_bv()));

  // Process loop
  AABB<S> local_aabb;
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task = task_stack.top();
    task_stack.pop();

    // Prune it if required
    if (try_prune_inner_nodes && (!this_task.is_leaf_node) &&
        prune_octree_info->prune_internal_nodes[this_task.node_vector_index]) {
      continue;
    }

    // If this is a leaf node
    if (this_task.is_leaf_node) {
      assert(this_task.node_vector_index < leaf_nodes.size());
      const OctreeLeafNode& leaf_node = leaf_nodes[this_task.node_vector_index];
      if (leaf_node.is_fully_occupied()) {
        const bool done = visitor(this_task.bv, this_task.depth, true);
        if (done) return;

        // Else, to the next task element
        continue;
      }

      // Into children
      assert(!leaf_node.is_fully_occupied());
      for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
        if (!leaf_node.child_occupied.test_i(child_i)) continue;
        computeChildAABB(this_task.bv, child_i, local_aabb);
        const bool done = visitor(local_aabb, this_task.depth + 1, true);
        if (done) return;
      }

      // End processing for leaf node
      continue;
    }

    // Not non-leaf case
    assert(!this_task.is_leaf_node);
    assert(this_task.node_vector_index < inner_nodes.size());
    const OctreeInnerNode& node = inner_nodes[this_task.node_vector_index];
    const bool node_full = inner_nodes_full[this_task.node_vector_index];

    // If this is a full node
    if (node_full) {
      const bool done = visitor(this_task.bv, this_task.depth, true);
      if (done) {
        return;
      }

      // To next node
      continue;
    }

    // Into children
    assert(!node_full);
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      // Query the child index
      const auto child_vector_index = node.children[child_i];
      if (child_vector_index == kInvalidNodeIndex) continue;

      // Into local bv
      computeChildAABB(this_task.bv, child_i, local_aabb);

      // Visit the inner node
      const bool done = visitor(this_task.bv, this_task.depth, false);
      if (done) return;

      // Push into stack
      auto child_frame = tree.makeStackElementChild(this_task, local_aabb,
                                                    child_vector_index);
      task_stack.push(std::move(child_frame));
    }
  }
}

template <typename S>
void visitOctree(const Octree<S>& tree, const VisitOctreeNodeFunc<S>& visitor) {
  visitOctree<S>(tree, nullptr, visitor);
}

}  // namespace octree2
}  // namespace fcl
