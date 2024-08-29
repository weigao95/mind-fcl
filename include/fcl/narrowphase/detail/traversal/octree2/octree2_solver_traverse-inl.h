//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
template <typename Shape>
void CollisionSolverOctree2<S>::octreeShapeIntersectImpl(
    const Octree2CollisionGeometry<S>& octree_geom, const Shape& shape,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_shape,
    OctreeLeafComputeCache& cache) const {
  // Gather the info
  if (octree_geom.raw_octree() == nullptr) return;
  const auto& inner_nodes = octree_geom.inner_nodes();
  const auto& inner_nodes_full = octree_geom.inner_nodes_fully_occupied();
  const auto& leaf_nodes = octree_geom.leaf_nodes();
  const auto* prune_internal_nodes = octree_geom.prune_internal_nodes();
  const bool try_prune_inner_nodes = (prune_internal_nodes != nullptr);
  assert(inner_nodes.size() == inner_nodes_full.size());

  // Make disjoint
  FixedRotationBoxDisjoint<S> disjoint;
  AABB<S> shape_local_AABB;
  {
    // Compute the bv for shape
    OBB<S> shape_obb_world;
    computeBV(shape, tf_shape, shape_obb_world);

    // Convert OBB to local AABB and a tf_AABB on that AABB
    shape_local_AABB.max_ = shape_obb_world.extent;
    shape_local_AABB.min_ = -shape_obb_world.extent;

    // Make tf_AABB frame
    Transform3<S> tf_shape_AABB;
    tf_shape_AABB.setIdentity();
    tf_shape_AABB.linear().matrix() = shape_obb_world.axis;
    tf_shape_AABB.translation() = shape_obb_world.To;
    disjoint.initialize(tf_octree, tf_shape_AABB);
  }

  // Make the stack
  using StackElement = octree2::OctreeTraverseStackElement<S>;
  std::stack<StackElement> task_stack;
  task_stack.push(StackElement::MakeRoot(octree_geom.octree_root_bv()));

  // Process loop
  AABB<S> local_aabb;
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task = task_stack.top();
    task_stack.pop();

    // Prune it if required
    if (try_prune_inner_nodes && (!this_task.is_leaf_node) &&
        prune_internal_nodes->operator[](this_task.node_vector_index)) {
      continue;
    }

    // Prune by obb
    const bool is_disjoint =
        disjoint.isDisjoint(this_task.bv, shape_local_AABB, false);
    if (is_disjoint) {
      continue;
    }

    // If this is a leaf
    if (this_task.is_leaf_node) {
      assert(this_task.node_vector_index < leaf_nodes.size());
      const auto& leaf_node = leaf_nodes[this_task.node_vector_index];
      if (leaf_node.is_fully_occupied()) {
        const auto encoded_node_idx =
            encodeOctree2Node(this_task.node_vector_index, true);
        boxToShapeProcessLeafPair<Shape>(octree_geom, tf_octree, shape,
                                         tf_shape, encoded_node_idx,
                                         this_task.bv, cache);
        if (request->terminationConditionSatisfied(*result)) return;
      } else {
        // Not fully occupied
        assert(!leaf_node.is_fully_occupied());
        for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
          if (!leaf_node.child_occupied.test_i(child_i)) continue;
          octree2::computeChildAABB(this_task.bv, child_i, local_aabb);

          // Invoke processor
          const auto encoded_node_idx =
              encodeOctree2Node(this_task.node_vector_index, true, child_i);
          boxToShapeProcessLeafPair<Shape>(octree_geom, tf_octree, shape,
                                           tf_shape, encoded_node_idx,
                                           local_aabb, cache);
          if (request->terminationConditionSatisfied(*result)) return;
        }
      }

      // Not finished
      continue;
    }

    // This is not leaf
    assert(!this_task.is_leaf_node);
    const octree2::OctreeInnerNode& node =
        inner_nodes[this_task.node_vector_index];
    const bool node_full = inner_nodes_full[this_task.node_vector_index];

    // If this is a full node
    if (node_full) {
      // Handle the node here
      const auto encoded_node_idx =
          encodeOctree2Node(this_task.node_vector_index, false);
      boxToShapeProcessLeafPair<Shape>(octree_geom, tf_octree, shape, tf_shape,
                                       encoded_node_idx, this_task.bv, cache);
      if (request->terminationConditionSatisfied(*result)) return;

      // To next node
      continue;
    }

    // Into children
    assert(!node_full);
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      const auto child_vector_index = node.children[child_i];
      if (child_vector_index == octree2::kInvalidNodeIndex) continue;

      // Into local bv amd push into stack
      octree2::computeChildAABB(this_task.bv, child_i, local_aabb);
      auto child_frame = octree_geom.makeStackElementChild(
          this_task, local_aabb, child_vector_index);
      task_stack.push(std::move(child_frame));
    }
  }
}

template <typename S>
template <typename BV>
void CollisionSolverOctree2<S>::octreeBVHIntersect(
    const Octree2CollisionGeometry<S>& octree_geom, const BVHModel<BV>* bvh,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_bvh,
    OctreeLeafComputeCache& cache) const {
  // Gather the info
  if (octree_geom.raw_octree() == nullptr) return;
  const auto& inner_nodes = octree_geom.inner_nodes();
  const auto& inner_nodes_full = octree_geom.inner_nodes_fully_occupied();
  const auto& leaf_nodes = octree_geom.leaf_nodes();
  const auto* prune_internal_nodes = octree_geom.prune_internal_nodes();
  const bool try_prune_inner_nodes = (prune_internal_nodes != nullptr);

  // The stack
  using OctreeNodeInStack = octree2::OctreeTraverseStackElement<S>;
  struct TaskStackElement {
    OctreeNodeInStack octree_node;
    int bvh_node;
  };

  // Make the root
  std::stack<TaskStackElement> task_stack;
  {
    TaskStackElement element;
    element.octree_node =
        OctreeNodeInStack::MakeRoot(octree_geom.octree_root_bv());
    element.bvh_node = 0;
    task_stack.push(std::move(element));
  }

  // Stack loop
  OBB<S> obb_1_octree;
  OBB<S> obb_2_bvh;
  AABB<S> local_aabb_octree;
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task = task_stack.top();
    task_stack.pop();
    const OctreeNodeInStack& octree_node_1 = this_task.octree_node;
    const int bvh_node_2 = this_task.bvh_node;

    // Prune octree node it if required
    if (try_prune_inner_nodes && (!octree_node_1.is_leaf_node) &&
        prune_internal_nodes->operator[](octree_node_1.node_vector_index)) {
      continue;
    }

    // Compute OBB
    const BV node_2_bv = bvh->getBV(bvh_node_2).bv;
    convertBV(octree_node_1.bv, tf_octree, obb_1_octree);
    convertBV(node_2_bv, tf_bvh, obb_2_bvh);
    if (!obb_1_octree.overlap(obb_2_bvh)) {
      continue;
    }

    // The correctness of this expression is subtle, inner_nodes_full can only
    // be accessed with inner node index, but when octree_node_1.is_leaf_node is
    // true, the later expression would not be evaluated
    const bool is_octree_node1_traverse_leaf =
        octree_node_1.is_leaf_node ||
        inner_nodes_full[octree_node_1.node_vector_index];

    // Leaf case
    if (bvh->getBV(bvh_node_2).isLeaf() && is_octree_node1_traverse_leaf) {
      // Obb overlap checked. Just obtain the triangle is supported
      const int primitive_id = bvh->getBV(bvh_node_2).primitiveId();
      const Simplex<S> simplex = bvh->getSimplex(primitive_id);

      // Not octree leaf, but a full node
      if (!octree_node_1.is_leaf_node) {
        assert(inner_nodes_full[octree_node_1.node_vector_index]);
        const auto encoded_node_idx =
            encodeOctree2Node(octree_node_1.node_vector_index, false);
        boxToSimplexProcessLeafPair(octree_geom, tf_octree, octree_node_1.bv,
                                    encoded_node_idx, bvh, tf_bvh, simplex,
                                    primitive_id, cache);
        if (request->terminationConditionSatisfied(*result)) return;
        continue;
      }

      // Try with octree leaf
      assert(octree_node_1.is_leaf_node);
      const octree2::OctreeLeafNode& leaf_node =
          leaf_nodes[octree_node_1.node_vector_index];
      if (leaf_node.is_fully_occupied()) {
        const auto encoded_node_idx =
            encodeOctree2Node(octree_node_1.node_vector_index, true);
        boxToSimplexProcessLeafPair(octree_geom, tf_octree, octree_node_1.bv,
                                    encoded_node_idx, bvh, tf_bvh, simplex,
                                    primitive_id, cache);
        if (request->terminationConditionSatisfied(*result)) return;
      } else {
        // Not fully occupied
        assert(!leaf_node.is_fully_occupied());
        for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
          if (!leaf_node.child_occupied.test_i(child_i)) continue;
          octree2::computeChildAABB(this_task.octree_node.bv, child_i,
                                    local_aabb_octree);

          const auto encoded_node_idx =
              encodeOctree2Node(octree_node_1.node_vector_index, true, child_i);
          boxToSimplexProcessLeafPair(octree_geom, tf_octree, local_aabb_octree,
                                      encoded_node_idx, bvh, tf_bvh, simplex,
                                      primitive_id, cache);
          if (request->terminationConditionSatisfied(*result)) return;
        }
      }

      // Finish with leaf case
      continue;
    }

    // Non-leaf case
    bool continue_on_bvh2 = is_octree_node1_traverse_leaf;
    continue_on_bvh2 |= ((!bvh->getBV(bvh_node_2).isLeaf()) &&
                         octree_node_1.bv.size() < node_2_bv.size());

    // Continue on bvh
    if (continue_on_bvh2) {
      TaskStackElement new_task_left;
      new_task_left.octree_node = this_task.octree_node;
      new_task_left.bvh_node = bvh->getBV(bvh_node_2).leftChild();
      task_stack.push(std::move(new_task_left));

      TaskStackElement new_task_right;
      new_task_right.octree_node = this_task.octree_node;
      new_task_right.bvh_node = bvh->getBV(bvh_node_2).rightChild();
      task_stack.push(std::move(new_task_right));
      continue;
    }

    // Continue on octree
    assert(!is_octree_node1_traverse_leaf);
    assert(!octree_node_1.is_leaf_node);
    assert(!inner_nodes_full[octree_node_1.node_vector_index]);
    const auto& inner_node = inner_nodes[octree_node_1.node_vector_index];
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      const auto child_vector_index = inner_node.children[child_i];
      if (child_vector_index == octree2::kInvalidNodeIndex) continue;

      // Push into stack
      octree2::computeChildAABB(octree_node_1.bv, child_i, local_aabb_octree);
      auto child_frame_octree = octree_geom.makeStackElementChild(
          octree_node_1, local_aabb_octree, child_vector_index);
      TaskStackElement new_stack_element;
      new_stack_element.octree_node = child_frame_octree;
      new_stack_element.bvh_node = this_task.bvh_node;
      task_stack.push(std::move(new_stack_element));
    }
  }
}

template <typename S>
void CollisionSolverOctree2<S>::octreePairIntersect(
    const Octree2CollisionGeometry<S>& tree1,
    const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
    const Transform3<S>& tf_tree2, const FixedRotationBoxDisjoint<S>& disjoint,
    OctreeLeafComputeCache& cache) const {
  // Gather the info
  if (tree1.raw_octree() == nullptr || tree2.raw_octree() == nullptr) return;
  const auto& inner_nodes_1 = tree1.inner_nodes();
  const auto& inner_nodes_full_1 = tree1.inner_nodes_fully_occupied();
  const auto& leaf_nodes_1 = tree1.leaf_nodes();
  const auto* prune_internal_nodes_1 = tree1.prune_internal_nodes();
  const bool try_prune_inner_nodes_1 = (prune_internal_nodes_1 != nullptr);

  const auto& inner_nodes_2 = tree2.inner_nodes();
  const auto& inner_nodes_full_2 = tree2.inner_nodes_fully_occupied();
  const auto& leaf_nodes_2 = tree2.leaf_nodes();
  const auto* prune_internal_nodes_2 = tree2.prune_internal_nodes();
  const bool try_prune_inner_nodes_2 = (prune_internal_nodes_2 != nullptr);

  // The stack
  using OctreeNodeInStack = octree2::OctreeTraverseStackElement<S>;
  struct TaskStackElement {
    OctreeNodeInStack tree1_node;
    OctreeNodeInStack tree2_node;
  };

  // Make the root
  std::stack<TaskStackElement> task_stack;
  {
    TaskStackElement element;
    element.tree1_node = OctreeNodeInStack::MakeRoot(tree1.octree_root_bv());
    element.tree2_node = OctreeNodeInStack::MakeRoot(tree2.octree_root_bv());
    task_stack.push(std::move(element));
  }

  // Stack loop
  AABB<S> local_aabb;
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskStackElement this_task = task_stack.top();
    task_stack.pop();
    const OctreeNodeInStack& tree1_node = this_task.tree1_node;
    const OctreeNodeInStack& tree2_node = this_task.tree2_node;

    // Prune octree node it if required
    if (try_prune_inner_nodes_1 && (!tree1_node.is_leaf_node) &&
        prune_internal_nodes_1->operator[](tree1_node.node_vector_index)) {
      continue;
    }

    if (try_prune_inner_nodes_2 && (!tree2_node.is_leaf_node) &&
        prune_internal_nodes_2->operator[](tree2_node.node_vector_index)) {
      continue;
    }

    // Check leaf and non-leaf case
    const bool is_tree1_traverse_leaf =
        tree1_node.is_leaf_node ||
        inner_nodes_full_1[tree1_node.node_vector_index];
    const bool is_tree2_traverse_leaf =
        tree2_node.is_leaf_node ||
        inner_nodes_full_2[tree2_node.node_vector_index];
    if (is_tree1_traverse_leaf && is_tree2_traverse_leaf) {
      // Both are leaf nodes
      if (tree1_node.is_leaf_node && tree2_node.is_leaf_node) {
        const auto& leaf_1 = leaf_nodes_1[tree1_node.node_vector_index];
        const auto& leaf_2 = leaf_nodes_2[tree2_node.node_vector_index];
        const bool finished =
            octreePairTwoLeafNode(tree1, tree2, tf_tree1, tf_tree2, tree1_node,
                                  tree2_node, leaf_1, leaf_2, disjoint, cache);
        if (finished) return;
        continue;
      }

      // Node2 is leaf node
      if ((!tree1_node.is_leaf_node) && tree2_node.is_leaf_node) {
        const auto& leaf_2 = leaf_nodes_2[tree2_node.node_vector_index];
        const bool finished = octreePairInnerNodeWithLeafNode(
            tree1, tree2, tf_tree1, tf_tree2, tree1_node, tree2_node, leaf_2,
            false, disjoint, cache);
        if (finished) return;
        continue;
      }

      // Node1 is leaf node
      if (tree1_node.is_leaf_node && (!tree2_node.is_leaf_node)) {
        const auto& leaf_1 = leaf_nodes_1[tree1_node.node_vector_index];
        const bool finished = octreePairInnerNodeWithLeafNode(
            tree2, tree1, tf_tree2, tf_tree1, tree2_node, tree1_node, leaf_1,
            true, disjoint, cache);
        if (finished) return;
        continue;
      }

      // Both are not leaf
      assert((!tree1_node.is_leaf_node) && (!tree2_node.is_leaf_node));
      octreePairInnerNodePairAsLeaf(tree1, tree2, tf_tree1, tf_tree2,
                                    tree1_node, tree2_node, disjoint, cache);
      if (request->terminationConditionSatisfied(*result)) return;
      continue;
    }

    // Not leaf, prune with obb check
    constexpr bool strict_obb_checking = false;
    if (disjoint.isDisjoint(tree1_node.bv, tree2_node.bv,
                            strict_obb_checking)) {
      continue;
    }

    // Branching strategy
    bool continue_on_tree2 = is_tree1_traverse_leaf;
    continue_on_tree2 |= ((!is_tree2_traverse_leaf) &&
                          tree1_node.bv.size() < tree2_node.bv.size());
    // Into child
    if (continue_on_tree2) {
      assert(!is_tree2_traverse_leaf);
      assert(!tree2_node.is_leaf_node);
      assert(!inner_nodes_full_2[tree2_node.node_vector_index]);
      const auto& inner_node = inner_nodes_2[tree2_node.node_vector_index];
      for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
        const auto child_vector_index = inner_node.children[child_i];
        if (child_vector_index == octree2::kInvalidNodeIndex) continue;

        // Push into stack
        octree2::computeChildAABB(tree2_node.bv, child_i, local_aabb);
        auto child_frame_octree = tree2.makeStackElementChild(
            tree2_node, local_aabb, child_vector_index);
        TaskStackElement new_stack_element;
        new_stack_element.tree1_node = tree1_node;
        new_stack_element.tree2_node = child_frame_octree;
        task_stack.push(std::move(new_stack_element));
      }

      // Done with node 2
      continue;
    }

    // Continue on node 1
    assert(!continue_on_tree2);
    assert(!is_tree1_traverse_leaf);
    assert(!tree1_node.is_leaf_node);
    assert(!inner_nodes_full_1[tree1_node.node_vector_index]);
    const auto& inner_node = inner_nodes_1[tree1_node.node_vector_index];
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      const auto child_vector_index = inner_node.children[child_i];
      if (child_vector_index == octree2::kInvalidNodeIndex) continue;

      // Push into stack
      octree2::computeChildAABB(tree1_node.bv, child_i, local_aabb);
      auto child_frame_octree = tree1.makeStackElementChild(
          tree1_node, local_aabb, child_vector_index);
      TaskStackElement new_stack_element;
      new_stack_element.tree1_node = child_frame_octree;
      new_stack_element.tree2_node = tree2_node;
      task_stack.push(std::move(new_stack_element));
    }
  }
}

}  // namespace detail
}  // namespace fcl
