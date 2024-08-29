#pragma once

namespace fcl {
namespace detail {

template <typename S>
template <typename Shape>
void TranslationalDisplacementOctreeSolver<S>::RunShapeOctree(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) const {
  // Check input
  if (s1 == nullptr || tree2 == nullptr || tree2->raw_octree() == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into impl
  runShapeOctree(s1, tf1, s1_displacement, tree2, tf2);
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementOctreeSolver<S>::RunOctreeShape(
    const Octree2CollisionGeometry<S>* tree1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& tree1_displacement, const Shape* s2,
    const Transform3<S>& tf2, const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) const {
  // Check input
  if (tree1 == nullptr || tree1->raw_octree() == nullptr || s2 == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into impl
  TranslationalDisplacement<S> displacement_tf;
  displacement_tf.scalar_displacement = tree1_displacement.scalar_displacement;
  displacement_tf.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf1.linear() * (-tree1_displacement.unit_axis_in_shape1));
  runShapeOctree(s2, tf2, displacement_tf, tree1, tf1);
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementOctreeSolver<S>::runShapeOctree(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2) const {
  // Gather the info
  const auto& octree_geom = *tree2;
  const auto& inner_nodes = octree_geom.inner_nodes();
  const auto& inner_nodes_full = octree_geom.inner_nodes_fully_occupied();
  const auto& leaf_nodes = octree_geom.leaf_nodes();
  const auto* prune_internal_nodes = octree_geom.prune_internal_nodes();
  const bool try_prune_inner_nodes = (prune_internal_nodes != nullptr);
  assert(inner_nodes.size() == inner_nodes_full.size());

  // Init the shape AABB and disjoint
  FixedOrientationBoxPairTranslationalCCD<S> disjoint;
  AABB<S> shape_local_AABB;
  initializeShapeFixedOrientationBoxTranslationalCCD<Shape>(
      *s1, tf1, s1_displacement, tf2, disjoint, shape_local_AABB);

  // Make the stack
  struct StackElement {
    octree2::OctreeTraverseStackElement<S> tree_node;
    Interval<S> interval;
  };

  // Init for task stack
  std::stack<StackElement> task_stack;
  {
    StackElement element;
    element.tree_node = octree2::OctreeTraverseStackElement<S>::MakeRoot(
        octree_geom.octree_root_bv());
    element.interval.lower_bound = 0.0;
    element.interval.upper_bound = 1.0;
    task_stack.push(std::move(element));
  }

  // Process loop
  AABB<S> local_aabb;
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task_and_interval = task_stack.top();
    task_stack.pop();
    const auto& this_task = this_task_and_interval.tree_node;
    const auto& bv_parent_interval = this_task_and_interval.interval;

    // Prune it if required
    if (try_prune_inner_nodes && (!this_task.is_leaf_node) &&
        prune_internal_nodes->operator[](this_task.node_vector_index)) {
      continue;
    }

    // Prune by obb
    Interval<S> obb_interval;
    const bool is_disjoint =
        disjoint.IsDisjoint(shape_local_AABB, this_task.bv, bv_parent_interval,
                            obb_interval, request->zero_movement_tolerance);
    if (is_disjoint) {
      continue;
    }

    // If this is a leaf
    if (this_task.is_leaf_node) {
      assert(this_task.node_vector_index < leaf_nodes.size());
      const auto& leaf_node = leaf_nodes[this_task.node_vector_index];
      if (leaf_node.is_fully_occupied()) {
        const auto encoded_node_idx =
            CollisionSolverOctree2<S>::encodeOctree2Node(
                this_task.node_vector_index, true);
        shapeToBoxProcessLeafPair<Shape>(s1, tf1, s1_displacement, tree2, tf2,
                                         encoded_node_idx, this_task.bv);

        // Check termination
        if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
          return;
        }
      } else {
        // Not fully occupied
        assert(!leaf_node.is_fully_occupied());
        for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
          if (!leaf_node.child_occupied.test_i(child_i)) continue;
          octree2::computeChildAABB(this_task.bv, child_i, local_aabb);

          // Invoke processor
          const auto encoded_node_idx =
              CollisionSolverOctree2<S>::encodeOctree2Node(
                  this_task.node_vector_index, true, child_i);
          shapeToBoxProcessLeafPair<Shape>(s1, tf1, s1_displacement, tree2, tf2,
                                           encoded_node_idx, local_aabb);

          // Check termination
          if (result->TerminationConditionSatisfied(
                  request->num_max_contacts)) {
            return;
          }
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
          CollisionSolverOctree2<S>::encodeOctree2Node(
              this_task.node_vector_index, false);
      shapeToBoxProcessLeafPair<Shape>(s1, tf1, s1_displacement, tree2, tf2,
                                       encoded_node_idx, this_task.bv);
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return;
      }

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
      StackElement child_frame;
      child_frame.tree_node = octree_geom.makeStackElementChild(
          this_task, local_aabb, child_vector_index);
      child_frame.interval = obb_interval;
      task_stack.push(std::move(child_frame));
    }
  }
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementOctreeSolver<S>::shapeToBoxProcessLeafPair(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2,
    std::int64_t encoded_octree_node_idx, const AABB<S>& voxel_aabb) const {
  // Make contact
  ContinuousContactMeta<S> contact;
  contact.reset();
  contact.o1 = s1;
  contact.o2 = tree2;
  contact.b1 = ContinuousContactMeta<S>::kNone;
  contact.b2 = encoded_octree_node_idx;
  contact.o2_bv = voxel_aabb;

  // Make box
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(voxel_aabb, tf2, box, box_tf);
  box.computeLocalAABB();

  // Invoke caller
  using Solver = ShapePairTranslationalCollisionSolver<S>;
  Solver::template RunShapePair<Shape, Box<S>>(
      *s1, tf1, s1_displacement, box, box_tf, contact, *request, *result);
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::RunOctreeObbBVH(
    const Octree2CollisionGeometry<S>* tree1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& tree1_displacement,
    const BVHModel<OBB<S>>* tree2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  if (tree1 == nullptr || tree2 == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Invoke the interface
  runOctreeObbBVH(tree1, tf1, tree1_displacement, tree2, tf2);
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::RunObbBVH_Octree(
    const BVHModel<OBB<S>>* tree1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& tree1_displacement,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  if (tree1 == nullptr || tree2 == nullptr) return;
  request = &request_in;
  result = &result_in;

  TranslationalDisplacement<S> displacement;
  displacement.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf1.linear() * (-tree1_displacement.unit_axis_in_shape1));
  displacement.scalar_displacement = tree1_displacement.scalar_displacement;
  runOctreeObbBVH(tree2, tf2, displacement, tree1, tf1);
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::runOctreeObbBVH(
    const Octree2CollisionGeometry<S>* tree1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const BVHModel<OBB<S>>* bvh2, const Transform3<S>& tf2) const {
  // Gather the info
  assert(tree1 != nullptr);
  const auto& octree_geom = *tree1;
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
    Interval<S> interval;
  };

  // Make the root
  std::stack<TaskStackElement> task_stack;
  {
    TaskStackElement element;
    element.octree_node =
        OctreeNodeInStack::MakeRoot(octree_geom.octree_root_bv());
    element.bvh_node = 0;
    element.interval.lower_bound = 0.0;
    element.interval.upper_bound = 1.0;
    task_stack.push(std::move(element));
  }

  // Stack loop
  OBB<S> obb1_octree;
  OBB<S> obb2_bvh;
  AABB<S> local_aabb_octree;
  while (!task_stack.empty()) {
    // Pop the stack
    const auto this_task = task_stack.top();
    task_stack.pop();
    const OctreeNodeInStack& octree_node_1 = this_task.octree_node;
    const int bvh_node_2 = this_task.bvh_node;
    const auto& bv_interval_parent = this_task.interval;

    // Prune octree node it if required
    if (try_prune_inner_nodes && (!octree_node_1.is_leaf_node) &&
        prune_internal_nodes->operator[](octree_node_1.node_vector_index)) {
      continue;
    }

    // Compute OBB
    const OBB<S>& node_2_raw_bv = bvh2->getBV(bvh_node_2).bv;
    convertBV(octree_node_1.bv, tf1, obb1_octree);
    convertBV(node_2_raw_bv, tf2, obb2_bvh);
    Interval<S> obb_interval;
    if (BoxPairTranslationalCCD<S>::IsDisjoint(
            obb1_octree, s1_displacement, obb2_bvh, bv_interval_parent,
            obb_interval, request->zero_movement_tolerance)) {
      continue;
    }

    // The correctness of this expression is subtle, inner_nodes_full can only
    // be accessed with inner node index, but when octree_node_1.is_leaf_node is
    // true, the later expression would not be evaluated
    const bool is_octree_node1_traverse_leaf =
        octree_node_1.is_leaf_node ||
        inner_nodes_full[octree_node_1.node_vector_index];

    // Leaf case
    if (bvh2->getBV(bvh_node_2).isLeaf() && is_octree_node1_traverse_leaf) {
      // Obb overlap checked. Just obtain the triangle is supported
      const int primitive_id = bvh2->getBV(bvh_node_2).primitiveId();
      const Simplex<S>& simplex = bvh2->getSimplex(primitive_id);

      // Not octree leaf, but a full node
      if (!octree_node_1.is_leaf_node) {
        assert(inner_nodes_full[octree_node_1.node_vector_index]);
        const auto encoded_node_idx =
            CollisionSolverOctree2<S>::encodeOctree2Node(
                octree_node_1.node_vector_index, false);
        boxToSimplexProcessLeafPair(tree1, tf1, s1_displacement,
                                    octree_node_1.bv, encoded_node_idx, bvh2,
                                    tf2, simplex, primitive_id);
        // Check termination
        if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
          return;
        }

        // Not terminated, continue
        continue;
      }

      // Try with octree leaf
      assert(octree_node_1.is_leaf_node);
      const octree2::OctreeLeafNode& leaf_node =
          leaf_nodes[octree_node_1.node_vector_index];
      if (leaf_node.is_fully_occupied()) {
        const auto encoded_node_idx =
            CollisionSolverOctree2<S>::encodeOctree2Node(
                octree_node_1.node_vector_index, true);
        boxToSimplexProcessLeafPair(tree1, tf1, s1_displacement,
                                    octree_node_1.bv, encoded_node_idx, bvh2,
                                    tf2, simplex, primitive_id);
        if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
          return;
        }
      } else {
        // Not fully occupied
        assert(!leaf_node.is_fully_occupied());
        for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
          if (!leaf_node.child_occupied.test_i(child_i)) continue;
          octree2::computeChildAABB(this_task.octree_node.bv, child_i,
                                    local_aabb_octree);

          const auto encoded_node_idx =
              CollisionSolverOctree2<S>::encodeOctree2Node(
                  octree_node_1.node_vector_index, true, child_i);
          boxToSimplexProcessLeafPair(tree1, tf1, s1_displacement,
                                      local_aabb_octree, encoded_node_idx, bvh2,
                                      tf2, simplex, primitive_id);
          if (result->TerminationConditionSatisfied(
                  request->num_max_contacts)) {
            return;
          }
        }
      }

      // Finish with leaf case
      continue;
    }

    // Non-leaf case
    bool continue_on_bvh2 = is_octree_node1_traverse_leaf;
    continue_on_bvh2 |= ((!bvh2->getBV(bvh_node_2).isLeaf()) &&
                         octree_node_1.bv.size() < node_2_raw_bv.size());

    // Continue on bvh
    if (continue_on_bvh2) {
      TaskStackElement new_task_left;
      new_task_left.octree_node = this_task.octree_node;
      new_task_left.bvh_node = bvh2->getBV(bvh_node_2).leftChild();
      new_task_left.interval = obb_interval;
      task_stack.push(std::move(new_task_left));

      TaskStackElement new_task_right;
      new_task_right.octree_node = this_task.octree_node;
      new_task_right.bvh_node = bvh2->getBV(bvh_node_2).rightChild();
      new_task_right.interval = obb_interval;
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
      new_stack_element.interval = obb_interval;
      task_stack.push(std::move(new_stack_element));
    }
  }
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::boxToSimplexProcessLeafPair(
    const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_octree,
    const TranslationalDisplacement<S>& octree_displacement,
    const AABB<S>& voxel_aabb, std::int64_t encoded_octree_node_idx,
    const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
    const Simplex<S>& simplex, fcl::intptr_t index_2) const {
  // Compute the box
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(voxel_aabb, tf_octree, box, box_tf);
  box.computeLocalAABB();

  // Make contact meta
  ContinuousContactMeta<S> contact;
  contact.o1 = octree;
  contact.o2 = bvh;
  contact.b1 = encoded_octree_node_idx;
  contact.b2 = index_2;
  contact.o1_bv = voxel_aabb;

  // Compute the contact
  using Solver = ShapePairTranslationalCollisionSolver<S>;
  Solver::template RunShapeSimplex<Box<S>>(box, box_tf, octree_displacement,
                                           simplex, tf_bvh, contact, *request,
                                           *result);
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::RunOctreePair(
    const Octree2CollisionGeometry<S>* tree1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& tree1_displacement,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  if (tree1 == nullptr || tree2 == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Init the shape AABB and disjoint
  FixedOrientationBoxPairTranslationalCCD<S> disjoint;
  disjoint.Initialize(tf1, tf2, tree1_displacement);
  runOctreePair(tree1, tree2, disjoint);
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::runOctreePair(
    const Octree2CollisionGeometry<S>* tree1_ptr,
    const Octree2CollisionGeometry<S>* tree2_ptr,
    const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const {
  // Gather the info
  assert(tree1_ptr != nullptr && tree2_ptr != nullptr);
  const auto& tree1 = *tree1_ptr;
  const auto& tree2 = *tree2_ptr;
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
    Interval<S> interval{};
  };

  // Make the root
  std::stack<TaskStackElement> task_stack;
  {
    TaskStackElement element;
    element.tree1_node = OctreeNodeInStack::MakeRoot(tree1.octree_root_bv());
    element.tree2_node = OctreeNodeInStack::MakeRoot(tree2.octree_root_bv());
    element.interval.lower_bound = 0.0;
    element.interval.upper_bound = 1.0;
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
    const Interval<S>& bv_parent_interval = this_task.interval;

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
            runOctreePairTwoLeafNode(tree1_ptr, tree2_ptr, tree1_node,
                                     tree2_node, leaf_1, leaf_2, disjoint);
        if (finished) return;
        continue;
      }

      // Node2 is leaf node
      if ((!tree1_node.is_leaf_node) && tree2_node.is_leaf_node) {
        const auto& leaf_2 = leaf_nodes_2[tree2_node.node_vector_index];
        const bool finished = runOctreePairInnerNodeWithLeafNode(
            tree1_ptr, tree2_ptr, tree1_node, tree2_node, leaf_2, false,
            disjoint);
        if (finished) return;
        continue;
      }

      // Node1 is leaf node
      if (tree1_node.is_leaf_node && (!tree2_node.is_leaf_node)) {
        const auto& leaf_1 = leaf_nodes_1[tree1_node.node_vector_index];
        const bool finished = runOctreePairInnerNodeWithLeafNode(
            tree2_ptr, tree1_ptr, tree2_node, tree1_node, leaf_1, true,
            disjoint);
        if (finished) return;
        continue;
      }

      // Both are not leaf
      assert((!tree1_node.is_leaf_node) && (!tree2_node.is_leaf_node));
      runOctreePairInnerNodePairAsLeaf(tree1_ptr, tree2_ptr, tree1_node,
                                       tree2_node, disjoint);
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return;
      }

      // Not terminated
      continue;
    }

    // Not leaf, prune with obb check
    Interval<S> obb_interval;
    if (disjoint.IsDisjoint(tree1_node.bv, tree2_node.bv, bv_parent_interval,
                            obb_interval, request->zero_movement_tolerance)) {
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
        new_stack_element.interval = obb_interval;
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
      new_stack_element.interval = obb_interval;
      task_stack.push(std::move(new_stack_element));
    }
  }
}

template <typename S>
void TranslationalDisplacementOctreeSolver<S>::runOctreePairInnerNodePairAsLeaf(
    const Octree2CollisionGeometry<S>* tree1,
    const Octree2CollisionGeometry<S>* tree2,
    const octree2::OctreeTraverseStackElement<S>& tree1_elem,
    const octree2::OctreeTraverseStackElement<S>& tree2_elem,
    const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const {
  ContinuousCollisionContact<S> contact;
  contact.o1 = tree1;
  contact.o2 = tree2;
  contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(
      tree1_elem.node_vector_index, false);
  contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
      tree2_elem.node_vector_index, false);
  contact.o1_bv = tree1_elem.bv;
  contact.o2_bv = tree2_elem.bv;

  // Invoke checking
  if (!disjoint.IsDisjoint(contact.o1_bv, contact.o2_bv, contact.toc,
                           request->zero_movement_tolerance)) {
    result->AddContact(contact);
  }
}

template <typename S>
bool TranslationalDisplacementOctreeSolver<S>::
    runOctreePairInnerNodeWithLeafNode(
        const Octree2CollisionGeometry<S>* tree1,
        const Octree2CollisionGeometry<S>* tree2,
        const octree2::OctreeTraverseStackElement<S>& tree1_elem,
        const octree2::OctreeTraverseStackElement<S>& tree2_elem,
        const octree2::OctreeLeafNode& tree2_node, bool reverse_tree12,
        const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const {
  ContinuousContactMeta<S> contact;
  contact.o1 = tree1;
  contact.o2 = tree2;
  contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(
      tree1_elem.node_vector_index, false);
  contact.o1_bv = tree1_elem.bv;
  contact.reverse_o1_and_o2 = reverse_tree12;

  // If this is a full node
  if (tree2_node.is_fully_occupied()) {
    // Assign the bv
    contact.o2_bv = tree2_elem.bv;
    contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
        tree2_elem.node_vector_index, true);

    // Check disjoint
    const bool is_disjoint =
        reverse_tree12 ? disjoint.IsDisjoint(contact.o2_bv, tree1_elem.bv,
                                             contact.external_box_toc,
                                             request->zero_movement_tolerance)
                       : disjoint.IsDisjoint(tree1_elem.bv, contact.o2_bv,
                                             contact.external_box_toc,
                                             request->zero_movement_tolerance);
    if (!is_disjoint) {
      ContinuousCollisionContact<S> actual_contact;
      contact.writeToContact(actual_contact);
      result->AddContact(std::move(actual_contact));
    }

    // Done
    return result->TerminationConditionSatisfied(request->num_max_contacts);
  }

  // Iterate over nodes in tree2
  assert(!tree2_node.is_fully_occupied());
  for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
    // Check occupied
    if (!tree2_node.child_occupied.test_i(child_i)) continue;

    // Always compute the child AABB
    octree2::computeChildAABB(tree2_elem.bv, child_i, contact.o2_bv);
    contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
        tree2_elem.node_vector_index, true, child_i);

    // Check disjoint
    const bool is_disjoint =
        reverse_tree12 ? disjoint.IsDisjoint(contact.o2_bv, tree1_elem.bv,
                                             contact.external_box_toc,
                                             request->zero_movement_tolerance)
                       : disjoint.IsDisjoint(tree1_elem.bv, contact.o2_bv,
                                             contact.external_box_toc,
                                             request->zero_movement_tolerance);
    if (!is_disjoint) {
      ContinuousCollisionContact<S> actual_contact;
      contact.writeToContact(actual_contact);
      result->AddContact(std::move(actual_contact));
    }

    // Done after update
    if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
      return true;
    }
  }

  // Not done
  return false;
}

template <typename S>
bool TranslationalDisplacementOctreeSolver<S>::runOctreePairTwoLeafNode(
    const Octree2CollisionGeometry<S>* tree1,
    const Octree2CollisionGeometry<S>* tree2,
    const octree2::OctreeTraverseStackElement<S>& tree1_elem,
    const octree2::OctreeTraverseStackElement<S>& tree2_elem,
    const octree2::OctreeLeafNode& tree1_node,
    const octree2::OctreeLeafNode& tree2_node,
    const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const {
  // Use full contact as cache
  ContinuousCollisionContact<S> contact;
  contact.o1 = tree1;
  contact.o2 = tree2;
  AABB<S>& o1_bv = contact.o1_bv;
  AABB<S>& o2_bv = contact.o2_bv;

  // Case 1: both are fully occupied
  if (tree1_node.is_fully_occupied() && tree2_node.is_fully_occupied()) {
    // clang-format off
    o1_bv = tree1_elem.bv;
    o2_bv = tree2_elem.bv;
    contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(tree1_elem.node_vector_index, true);
    contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(tree2_elem.node_vector_index, true);
    // clang-format on

    // Invoke checking
    if (!disjoint.IsDisjoint(o1_bv, o2_bv, contact.toc,
                             request->zero_movement_tolerance)) {
      result->AddContact(std::move(contact));
    }

    // Done
    return result->TerminationConditionSatisfied(request->num_max_contacts);
  }

  // Case 2:
  if (tree1_node.is_fully_occupied() && (!tree2_node.is_fully_occupied())) {
    // Assign the bv for o1
    o1_bv = tree1_elem.bv;
    contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(
        tree1_elem.node_vector_index, true);

    // For another child
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      if (!tree2_node.child_occupied.test_i(child_i)) continue;
      octree2::computeChildAABB(tree2_elem.bv, child_i, o2_bv);
      contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
          tree2_elem.node_vector_index, true, child_i);

      // Invoke checking
      if (!disjoint.IsDisjoint(o1_bv, o2_bv, contact.toc,
                               request->zero_movement_tolerance)) {
        result->AddContact(contact);  // Cannot move as it is used as cache
      }

      // Done with this case
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return true;
      }
    }

    // Not completely finished
    return false;
  }

  // Case 3
  if (tree2_node.is_fully_occupied() && (!tree1_node.is_fully_occupied())) {
    // Assign the bv for o2
    o2_bv = tree2_elem.bv;
    contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
        tree2_elem.node_vector_index, true);

    // Iterate the children
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      if (!tree1_node.child_occupied.test_i(child_i)) continue;
      octree2::computeChildAABB(tree1_elem.bv, child_i, o1_bv);
      contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(
          tree1_elem.node_vector_index, true, child_i);

      // Invoke checking
      if (!disjoint.IsDisjoint(o1_bv, o2_bv, contact.toc,
                               request->zero_movement_tolerance)) {
        result->AddContact(contact);  // Cannot move as it is used as cache
      }

      // Done with this case
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return true;
      }
    }

    // Not completely finished
    return false;
  }

  // Both are not fully satisfied
  assert(!tree1_node.is_fully_occupied());
  assert(!tree2_node.is_fully_occupied());
  if (disjoint.IsDisjoint(tree1_elem.bv, tree2_elem.bv, contact.toc,
                          request->zero_movement_tolerance)) {
    return false;
  }

  // Might prune the children
  octree2::Bitset8 tree1_child2check = tree1_node.child_occupied;
  octree2::Bitset8 tree2_child2check = tree2_node.child_occupied;
  {
    const int tree1_n_child = tree1_child2check.count_of_bits();
    const int tree2_n_child = tree2_child2check.count_of_bits();
    const bool use_parent_obb_prune =
        tree1_n_child * tree2_n_child > tree1_n_child + tree2_n_child;
    if (use_parent_obb_prune) {
      // Prune tree1 child
      for (std::uint8_t i = 0; i < 8; i++) {
        if (!tree1_child2check.test_i(i)) continue;
        octree2::computeChildAABB(tree1_elem.bv, i, o1_bv);
        const bool is_disjoint =
            disjoint.IsDisjoint(o1_bv, tree2_elem.bv, contact.toc,
                                request->zero_movement_tolerance);
        if (is_disjoint) tree1_child2check.clear_i(i);
      }

      // Prune tree2 child
      for (std::uint8_t i = 0; i < 8; i++) {
        if (!tree2_child2check.test_i(i)) continue;
        octree2::computeChildAABB(tree2_elem.bv, i, o2_bv);
        const bool is_disjoint =
            disjoint.IsDisjoint(tree1_elem.bv, o2_bv, contact.toc,
                                request->zero_movement_tolerance);
        if (is_disjoint) tree2_child2check.clear_i(i);
      }
    }
  }

  // Loop over child
  for (std::uint8_t child_1 = 0; child_1 < 8; child_1++) {
    if (!tree1_child2check.test_i(child_1)) continue;
    octree2::computeChildAABB(tree1_elem.bv, child_1, o1_bv);
    contact.b1 = CollisionSolverOctree2<S>::encodeOctree2Node(
        tree1_elem.node_vector_index, true, child_1);

    // Into child 2
    for (std::uint8_t child_2 = 0; child_2 < 8; child_2++) {
      if (!tree2_child2check.test_i(child_2)) continue;
      octree2::computeChildAABB(tree2_elem.bv, child_2, o2_bv);
      contact.b2 = CollisionSolverOctree2<S>::encodeOctree2Node(
          tree2_elem.node_vector_index, true, child_2);

      // Invoke checking
      if (!disjoint.IsDisjoint(o1_bv, o2_bv, contact.toc,
                               request->zero_movement_tolerance)) {
        result->AddContact(contact);  // Cannot move as it is used as cache
      }

      // Done after update
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return true;
      }
    }
  }

  // Done
  return false;
}

}  // namespace detail
}  // namespace fcl