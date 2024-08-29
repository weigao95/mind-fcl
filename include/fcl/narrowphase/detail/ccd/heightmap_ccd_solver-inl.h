#pragma once

namespace fcl {
namespace detail {

// Shape with HeightMap
//==============================================================================
template <typename S>
template <typename Shape>
void TranslationalDisplacementHeightMapSolver<S>::runShapeHeightMap(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2) const {
  using heightmap::Pixel;
  const LayeredHeightMap& heightmap = *(hm2->raw_heightmap());

  // Init the shape AABB and disjoint
  FixedOrientationBoxPairTranslationalCCD<S> disjoint;
  AABB<S> shape_local_AABB;
  initializeShapeFixedOrientationBoxTranslationalCCD<Shape>(
      *s1, tf1, s1_displacement, tf2, disjoint, shape_local_AABB);

  // Define the task frame
  struct TaskFrame {
    HeightMapNode hm_node{};
    Interval<S> interval{};

    explicit TaskFrame() = default;
    explicit TaskFrame(HeightMapNode hm_node_in, Interval<S> interval_in)
        : hm_node(std::move(hm_node_in)), interval(std::move(interval_in)) {}
    ~TaskFrame() = default;
  };
  std::stack<TaskFrame> task_stack;

  // Init the task stack
  {
    const FlatHeightMap& hm_top = heightmap.top();
    for (uint16_t y = 0; y < hm_top.full_shape_y(); y++) {
      for (uint16_t x = 0; x < hm_top.full_shape_x(); x++) {
        TaskFrame task_xy;
        task_xy.hm_node = HeightMapNode(0, Pixel(x, y));
        task_xy.interval.lower_bound = 0.0;
        task_xy.interval.upper_bound = 1.0;
        task_stack.push(std::move(task_xy));
      }
    }
  }

  // The processing loop
  AABB<S> hm_aabb_local;
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskFrame this_task = task_stack.top();
    const HeightMapNode& hm_node = this_task.hm_node;
    const auto& bv_interval_parent = this_task.interval;
    task_stack.pop();

    // Obtain heightmap aabb
    if (!heightMapNodeToAABB<S>(heightmap, hm_node, hm_aabb_local)) {
      continue;
    }

    // OBB check
    Interval<S> obb_interval;
    if (disjoint.IsDisjoint(shape_local_AABB, hm_aabb_local, bv_interval_parent,
                            obb_interval, request->zero_movement_tolerance)) {
      continue;
    }

    if (heightmap.is_bottom_layer(hm_node.layer)) {
      // Both are leaf nodes
      shapeToBoxProcessLeafPair(s1, tf1, s1_displacement, hm2, tf2,
                                hm_node.pixel, hm_aabb_local);
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return;
      }
    } else {
      // clang-format off
      task_stack.push(TaskFrame(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 0, hm_node.pixel.y * 2 + 0)}, obb_interval));
      task_stack.push(TaskFrame(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 1, hm_node.pixel.y * 2 + 0)}, obb_interval));
      task_stack.push(TaskFrame(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 0, hm_node.pixel.y * 2 + 1)}, obb_interval));
      task_stack.push(TaskFrame(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 1, hm_node.pixel.y * 2 + 1)}, obb_interval));
      // clang-format on
    }
  }
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementHeightMapSolver<S>::shapeToBoxProcessLeafPair(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
    const heightmap::Pixel& pixel, const AABB<S>& pixel_aabb) const {
  ContinuousContactMeta<S> contact;
  contact.o1 = s1;
  contact.b1 = ContinuousContactMeta<S>::kNone;
  contact.o2 = hm2;
  contact.b2 = encodePixel(pixel);
  contact.o2_bv = pixel_aabb;

  // Make box
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(pixel_aabb, tf2, box, box_tf);
  box.computeLocalAABB();

  // Invoke caller
  using Solver = ShapePairTranslationalCollisionSolver<S>;
  Solver::template RunShapePair<Shape, Box<S>>(
      *s1, tf1, s1_displacement, box, box_tf, contact, *request, *result);
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementHeightMapSolver<S>::RunShapeHeightMap(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (s1 == nullptr || hm2 == nullptr || hm2->raw_heightmap() == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into impl
  runShapeHeightMap(s1, tf1, s1_displacement, hm2, tf2);
}

template <typename S>
template <typename Shape>
void TranslationalDisplacementHeightMapSolver<S>::RunHeightMapShape(
    const HeightMapCollisionGeometry<S>* hm1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement, const Shape* s2,
    const Transform3<S>& tf2, const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (s2 == nullptr || hm1 == nullptr || hm1->raw_heightmap() == nullptr) {
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
  displacement_tf.scalar_displacement = hm1_displacement.scalar_displacement;
  displacement_tf.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf1.linear() * (-hm1_displacement.unit_axis_in_shape1));
  runShapeHeightMap(s2, tf2, displacement_tf, hm1, tf1);
}

// BVH with HeightMap
//==============================================================================
template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::runHeightMapObbBVH(
    const HeightMapCollisionGeometry<S>* hm_geometry,
    const Transform3<S>& tf_hm,
    const TranslationalDisplacement<S>& hm1_displacement,
    const BVHModel<OBB<S>>* bvh, const Transform3<S>& tf_bvh) const {
  assert(hm_geometry != nullptr);
  assert(bvh != nullptr);

  // We need a stack following the recursive procedure
  struct TaskFrame {
    HeightMapNode hm_node;
    int bvh_node;
    Interval<S> interval;
  };

  // Init the task stack
  using heightmap::Pixel;
  std::stack<TaskFrame> task_stack;
  const heightmap::LayeredHeightMap<S>& hm_1 = *hm_geometry->raw_heightmap();
  {
    const FlatHeightMap& hm_top = hm_1.top();
    for (uint16_t hm1_y = 0; hm1_y < hm_top.full_shape_y(); hm1_y++) {
      for (uint16_t hm1_x = 0; hm1_x < hm_top.full_shape_x(); hm1_x++) {
        TaskFrame this_task;
        this_task.hm_node = HeightMapNode(0, Pixel(hm1_x, hm1_y));
        this_task.bvh_node = 0;
        this_task.interval.lower_bound = 0.0;
        this_task.interval.upper_bound = 1.0;
        task_stack.push(std::move(this_task));
      }
    }
  }

  // Stack loop
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskFrame this_task = task_stack.top();
    task_stack.pop();
    const HeightMapNode& hm_node_1 = this_task.hm_node;
    const int bvh_node_2 = this_task.bvh_node;
    const auto& bv_parent_interval = this_task.interval;

    // Compute the bv
    AABB<S> hm_aabb_1;
    const bool hm_aabb_1_valid =
        heightMapNodeToAABB<S>(hm_1, hm_node_1, hm_aabb_1);
    if (!hm_aabb_1_valid) continue;
    const OBB<S>& node_2_bv = bvh->getBV(bvh_node_2).bv;

    // Compute obb
    OBB<S> obb_1_hm, obb_2_bvh;
    convertBV(hm_aabb_1, tf_hm, obb_1_hm);
    convertBV(node_2_bv, tf_bvh, obb_2_bvh);
    Interval<S> obb_interval;
    if (BoxPairTranslationalCCD<S>::IsDisjoint(
            obb_1_hm, hm1_displacement, obb_2_bvh, bv_parent_interval,
            obb_interval, request->zero_movement_tolerance)) {
      continue;
    }

    // The root case
    if (bvh->getBV(bvh_node_2).isLeaf() &&
        hm_1.is_bottom_layer(hm_node_1.layer)) {
      // Obb overlap checked. Just obtain the triangle is supported
      const int primitive_id = bvh->getBV(bvh_node_2).primitiveId();
      const Simplex<S> simplex = bvh->getSimplex(primitive_id);
      boxToSimplexProcessLeafPair(hm_geometry, tf_hm, hm1_displacement,
                                  hm_aabb_1,
                                  heightmap::encodePixel(hm_node_1.pixel), bvh,
                                  tf_bvh, simplex, primitive_id);

      // Finished
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return;
      }
    } else {
      // Obb overlap checked. At least one of them is not the leaf
      bool continue_on_hm1 = bvh->getBV(bvh_node_2).isLeaf();
      continue_on_hm1 |= ((!hm_1.is_bottom_layer(hm_node_1.layer)) &&
                          hm_aabb_1.size() > node_2_bv.size());
      if (continue_on_hm1) {
        // Continue on hm_1
        assert(!hm_1.is_bottom_layer(hm_node_1.layer));
        std::array<HeightMapNode, 4> node_1_children;
        // clang-format off
        node_1_children[0] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 0, hm_node_1.pixel.y * 2 + 0)};
        node_1_children[1] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 1, hm_node_1.pixel.y * 2 + 0)};
        node_1_children[2] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 0, hm_node_1.pixel.y * 2 + 1)};
        node_1_children[3] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 1, hm_node_1.pixel.y * 2 + 1)};
        // clang-format on

        // Add the tasks
        for (std::size_t i = 0; i < 4; i++) {
          TaskFrame task_i;
          task_i.hm_node = node_1_children[i];
          task_i.bvh_node = bvh_node_2;
          task_i.interval = obb_interval;
          task_stack.push(std::move(task_i));
        }
      } else {
        // Continue on bvh_2, push left child to stack
        TaskFrame new_task_left;
        new_task_left.hm_node = hm_node_1;
        new_task_left.bvh_node = bvh->getBV(bvh_node_2).leftChild();
        new_task_left.interval = obb_interval;
        task_stack.push(std::move(new_task_left));

        // Push right child to stack
        TaskFrame new_task_right;
        new_task_right.hm_node = hm_node_1;
        new_task_right.bvh_node = bvh->getBV(bvh_node_2).rightChild();
        new_task_right.interval = obb_interval;
        task_stack.push(std::move(new_task_right));
      }
    }
  }
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::boxToSimplexProcessLeafPair(
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_hm,
    const TranslationalDisplacement<S>& hm1_displacement,
    const AABB<S>& node_aabb_hm, std::int64_t index_1,
    const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
    const Simplex<S>& simplex, std::int64_t index_2) const {
  // Compoute box
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(node_aabb_hm, tf_hm, box, box_tf);
  box.computeLocalAABB();

  // Make the meta and shape solver
  ContinuousContactMeta<S> contact;
  contact.o1 = heightmap_geometry;
  contact.o2 = bvh;
  contact.b1 = index_1;
  contact.b2 = index_2;
  contact.o1_bv = node_aabb_hm;

  // Compute the contact
  using Solver = ShapePairTranslationalCollisionSolver<S>;
  Solver::template RunShapeSimplex<Box<S>>(box, box_tf, hm1_displacement,
                                           simplex, tf_bvh, contact, *request,
                                           *result);
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::RunHeightMapObbBVH(
    const HeightMapCollisionGeometry<S>* hm1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement,
    const BVHModel<OBB<S>>* bvh2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (hm1 == nullptr || hm1->raw_heightmap() == nullptr || bvh2 == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into impl
  runHeightMapObbBVH(hm1, tf1, hm1_displacement, bvh2, tf2);
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::RunObbBVH_HeightMap(
    const BVHModel<OBB<S>>* bvh1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& bvh1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (hm2 == nullptr || hm2->raw_heightmap() == nullptr || bvh1 == nullptr) {
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
  displacement_tf.scalar_displacement = bvh1_displacement.scalar_displacement;
  displacement_tf.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf1.linear() * (-bvh1_displacement.unit_axis_in_shape1));
  runHeightMapObbBVH(hm2, tf2, displacement_tf, bvh1, tf1);
}

// Octree with HeightMap
//==============================================================================
template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::runHeightMapOctree(
    const HeightMapCollisionGeometry<S>* hm_geometry,
    const Transform3<S>& tf_hm,
    const TranslationalDisplacement<S>& hm_displacement,
    const Octree2CollisionGeometry<S>* octree,
    const Transform3<S>& tf_octree) const {
  assert(hm_geometry != nullptr);
  assert(octree != nullptr);
  if (hm_geometry == nullptr || octree == nullptr) return;

  // Obtain info for octree
  if (octree->raw_octree() == nullptr) return;
  const auto& inner_nodes = octree->inner_nodes();
  const auto& inner_nodes_full = octree->inner_nodes_fully_occupied();
  const auto& leaf_nodes = octree->leaf_nodes();
  const auto* prune_internal_nodes = octree->prune_internal_nodes();
  const bool try_prune_inner_nodes = (prune_internal_nodes != nullptr);

  // Define the stack
  using heightmap::Pixel;
  using OctreeNodeInStack = octree2::OctreeTraverseStackElement<S>;
  struct TaskStackElement {
    HeightMapNode hm_node_1;
    std::uint8_t depth_from_root_1{0};
    OctreeNodeInStack octree_node_2;
    Interval<S> interval;
  };

  // Init the task stack
  std::stack<TaskStackElement> task_stack;
  const LayeredHeightMap& hm_1 = *hm_geometry->raw_heightmap();
  {
    const FlatHeightMap& hm_top = hm_1.top();
    for (uint16_t hm1_y = 0; hm1_y < hm_top.full_shape_y(); hm1_y++) {
      for (uint16_t hm1_x = 0; hm1_x < hm_top.full_shape_x(); hm1_x++) {
        TaskStackElement this_task;
        this_task.hm_node_1 = HeightMapNode(0, Pixel(hm1_x, hm1_y));
        this_task.depth_from_root_1 = 0;
        this_task.octree_node_2 =
            OctreeNodeInStack::MakeRoot(octree->octree_root_bv());
        this_task.interval.lower_bound = 0.0;
        this_task.interval.upper_bound = 1.0;
        task_stack.push(std::move(this_task));
      }
    }
  }

  // Init the shape AABB and disjoint
  FixedOrientationBoxPairTranslationalCCD<S> disjoint;
  disjoint.Initialize(tf_hm, tf_octree, hm_displacement);

  // Task processing loop
  AABB<S> aabb_1_hm, child_aabb_tree_2;
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskStackElement this_task = task_stack.top();
    task_stack.pop();

    // Obtain the task data
    const auto& hm_node_1 = this_task.hm_node_1;
    const auto& octree_node_2 = this_task.octree_node_2;
    const auto& bv_interval_parent = this_task.interval;

    // Prune octree node it if required
    if (try_prune_inner_nodes && (!octree_node_2.is_leaf_node) &&
        prune_internal_nodes->operator[](octree_node_2.node_vector_index)) {
      continue;
    }

    // Check validity of hm
    const bool aabb_1_valid =
        heightMapNodeToAABB<S>(hm_1, hm_node_1, aabb_1_hm);
    if (!aabb_1_valid) continue;

    // Leaf condition of octree
    const bool is_octree_node2_traverse_leaf =
        octree_node_2.is_leaf_node ||
        inner_nodes_full[octree_node_2.node_vector_index];
    const bool is_hm_node1_bottom = hm_1.is_bottom_layer(hm_node_1.layer);
    if (is_octree_node2_traverse_leaf && is_hm_node1_bottom) {
      // This is inner node as leaf
      if (!octree_node_2.is_leaf_node) {
        assert(inner_nodes_full[octree_node_2.node_vector_index]);
        boxToBoxProcessLeafPair(hm_geometry, aabb_1_hm,
                                heightmap::encodePixel(hm_node_1.pixel), octree,
                                octree_node_2.bv,
                                octree_node_2.node_vector_index, disjoint);

        // Check termination
        if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
          return;
        }

        // Not terminated
        continue;
      }

      // This is real leaf
      assert(octree_node_2.is_leaf_node);
      const auto& leaf_2 = leaf_nodes[octree_node_2.node_vector_index];
      if (leaf_2.is_fully_occupied()) {
        boxToBoxProcessLeafPair(hm_geometry, aabb_1_hm,
                                heightmap::encodePixel(hm_node_1.pixel), octree,
                                octree_node_2.bv,
                                octree_node_2.node_vector_index, disjoint);

        // Check termination
        if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
          return;
        }
      } else {
        for (std::uint8_t i = 0; i < 8; i++) {
          if (!leaf_2.child_occupied.test_i(i)) continue;
          octree2::computeChildAABB(octree_node_2.bv, i, child_aabb_tree_2);
          boxToBoxProcessLeafPair(hm_geometry, aabb_1_hm,
                                  heightmap::encodePixel(hm_node_1.pixel),
                                  octree, child_aabb_tree_2,
                                  octree_node_2.node_vector_index, disjoint);

          // Check termination
          if (result->TerminationConditionSatisfied(
                  request->num_max_contacts)) {
            return;
          }
        }
      }

      // Finish the leaf case
      continue;
    }

    // Non-leaf case
    Interval<S> obb_interval;
    if (disjoint.IsDisjoint(aabb_1_hm, octree_node_2.bv, bv_interval_parent,
                            obb_interval, request->zero_movement_tolerance)) {
      continue;
    }

    // TODO(wei): inspect this parameter and overall branching strategy,
    //  which has a major impact on speed
    const std::uint8_t prefer_octree_offset = 0;
    bool continue_on_hm1 = is_octree_node2_traverse_leaf;
    continue_on_hm1 |= (!is_hm_node1_bottom) &&
                       (this_task.depth_from_root_1 + prefer_octree_offset <
                        this_task.octree_node_2.depth);
    if (continue_on_hm1) {
      assert(!is_hm_node1_bottom);
      std::array<HeightMapNode, 4> node_1_children;
      // clang-format off
      node_1_children[0] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 0, hm_node_1.pixel.y * 2 + 0)};
      node_1_children[1] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 1, hm_node_1.pixel.y * 2 + 0)};
      node_1_children[2] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 0, hm_node_1.pixel.y * 2 + 1)};
      node_1_children[3] = HeightMapNode{hm_node_1.layer + 1, Pixel(hm_node_1.pixel.x * 2 + 1, hm_node_1.pixel.y * 2 + 1)};
      // clang-format on

      // Add the tasks
      for (std::size_t i = 0; i < 4; i++) {
        TaskStackElement task_i;
        task_i.hm_node_1 = node_1_children[i];
        task_i.depth_from_root_1 = this_task.depth_from_root_1 + 1;
        task_i.octree_node_2 = this_task.octree_node_2;
        task_i.interval = obb_interval;
        task_stack.push(std::move(task_i));
      }

      // End for hm1 expansion
      continue;
    }

    // Not on hm1
    assert(!continue_on_hm1);
    assert(!is_octree_node2_traverse_leaf);
    assert(!octree_node_2.is_leaf_node);
    const auto& node = inner_nodes[octree_node_2.node_vector_index];
    for (std::uint8_t i = 0; i < 8; i++) {
      const std::uint32_t child_vector_index = node.children[i];
      if (child_vector_index == octree2::kInvalidNodeIndex) continue;

      // Push into stack
      octree2::computeChildAABB(octree_node_2.bv, i, child_aabb_tree_2);
      auto child_frame_octree = octree->makeStackElementChild(
          octree_node_2, child_aabb_tree_2, child_vector_index);
      TaskStackElement new_stack_element;
      new_stack_element.hm_node_1 = this_task.hm_node_1;
      new_stack_element.octree_node_2 = child_frame_octree;
      new_stack_element.interval = obb_interval;
      task_stack.push(std::move(new_stack_element));
    }
  }
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::boxToBoxProcessLeafPair(
    const fcl::CollisionGeometry<S>* geom_1, const AABB<S>& node_aabb_1,
    std::int64_t index_1, const fcl::CollisionGeometry<S>* geom_2,
    const AABB<S>& node_aabb_2, std::int64_t index_2,
    const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const {
  ContinuousCollisionContact<S> contact;
  contact.o1 = geom_1;
  contact.o2 = geom_2;
  contact.b1 = index_1;
  contact.b2 = index_2;
  contact.o1_bv = node_aabb_1;
  contact.o2_bv = node_aabb_2;

  if (!disjoint.IsDisjoint(contact.o1_bv, contact.o2_bv, contact.toc,
                           request->zero_movement_tolerance)) {
    result->AddContact(std::move(contact));
  }
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::RunHeightMapOctree(
    const HeightMapCollisionGeometry<S>* hm1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement,
    const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (hm1 == nullptr || hm1->raw_heightmap() == nullptr || octree == nullptr ||
      octree->raw_octree() == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into interface
  runHeightMapOctree(hm1, tf1, hm1_displacement, octree, tf2);
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::RunOctreeHeightMap(
    const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_octree,
    const TranslationalDisplacement<S>& octree_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (hm2 == nullptr || hm2->raw_heightmap() == nullptr || octree == nullptr ||
      octree->raw_octree() == nullptr) {
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
  displacement_tf.scalar_displacement = octree_displacement.scalar_displacement;
  displacement_tf.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf_octree.linear() * (-octree_displacement.unit_axis_in_shape1));
  runHeightMapOctree(hm2, tf2, displacement_tf, octree, tf_octree);
}

// HeightMap Pair
//==============================================================================
template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::runHeightMapPair(
    const HeightMapCollisionGeometry<S>* hm1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2) const {
  // The task stack contains the element for checking
  assert(hm1 != nullptr);
  assert(hm2 != nullptr);
  using heightmap::Pixel;
  const auto& hm_1_geometry = *hm1;
  const auto& hm_2_geometry = *hm2;
  const LayeredHeightMap& hm_1 = *hm_1_geometry.raw_heightmap();
  const LayeredHeightMap& hm_2 = *hm_2_geometry.raw_heightmap();
  struct TaskFrame {
    HeightMapNode node_1;
    HeightMapNode node_2;
    Interval<S> interval;
  };
  std::stack<TaskFrame> task_stack;

  // Init the task stack
  {
    const FlatHeightMap& hm1_top = hm_1.top();
    const FlatHeightMap& hm2_top = hm_2.top();
    for (uint16_t hm1_y = 0; hm1_y < hm1_top.full_shape_y(); hm1_y++) {
      for (uint16_t hm1_x = 0; hm1_x < hm1_top.full_shape_x(); hm1_x++) {
        for (uint16_t hm2_y = 0; hm2_y < hm2_top.full_shape_y(); hm2_y++) {
          for (uint16_t hm2_x = 0; hm2_x < hm2_top.full_shape_x(); hm2_x++) {
            TaskFrame this_task;
            this_task.node_1 = HeightMapNode(0, Pixel(hm1_x, hm1_y));
            this_task.node_2 = HeightMapNode(0, Pixel(hm2_x, hm2_y));
            this_task.interval.lower_bound = 0.0;
            this_task.interval.upper_bound = 1.0;
            task_stack.push(std::move(this_task));
          }
        }
      }
    }
  }

  // Init the shape AABB and disjoint
  FixedOrientationBoxPairTranslationalCCD<S> disjoint;
  disjoint.Initialize(tf1, tf2, hm1_displacement);

  // The processing loop
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskFrame this_task = task_stack.top();
    const HeightMapNode& node_1 = this_task.node_1;
    const HeightMapNode& node_2 = this_task.node_2;
    const auto& bv_interval_parent = this_task.interval;
    task_stack.pop();

    // Obtain AABB
    AABB<S> aabb_1, aabb_2;
    const bool aabb_1_valid = heightMapNodeToAABB<S>(hm_1, node_1, aabb_1);
    const bool aabb_2_valid = heightMapNodeToAABB<S>(hm_2, node_2, aabb_2);
    if ((!aabb_1_valid) || (!aabb_2_valid)) continue;

    if (hm_1.is_bottom_layer(node_1.layer) &&
        hm_2.is_bottom_layer(node_2.layer)) {
      // Both are leaf nodes
      boxToBoxProcessLeafPair(hm1, aabb_1, heightmap::encodePixel(node_1.pixel),
                              hm2, aabb_2, heightmap::encodePixel(node_2.pixel),
                              disjoint);
      if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
        return;
      }
    } else {
      Interval<S> obb_interval;
      if (disjoint.IsDisjoint(aabb_1, aabb_2, bv_interval_parent, obb_interval,
                              request->zero_movement_tolerance)) {
        continue;
      }

      // Need to continue on box
      bool continue_on_1 = hm_2.is_bottom_layer(node_2.layer);
      continue_on_1 |= ((!hm_1.is_bottom_layer(node_1.layer)) &&
                        node_1.layer < node_2.layer);
      if (continue_on_1) {
        // Continue on hm_1
        assert(!hm_1.is_bottom_layer(node_1.layer));
        std::array<HeightMapNode, 4> node_1_children;
        // clang-format off
        node_1_children[0] = HeightMapNode{node_1.layer + 1, Pixel(node_1.pixel.x * 2 + 0, node_1.pixel.y * 2 + 0)};
        node_1_children[1] = HeightMapNode{node_1.layer + 1, Pixel(node_1.pixel.x * 2 + 1, node_1.pixel.y * 2 + 0)};
        node_1_children[2] = HeightMapNode{node_1.layer + 1, Pixel(node_1.pixel.x * 2 + 0, node_1.pixel.y * 2 + 1)};
        node_1_children[3] = HeightMapNode{node_1.layer + 1, Pixel(node_1.pixel.x * 2 + 1, node_1.pixel.y * 2 + 1)};
        // clang-format on

        // Add the tasks
        for (std::size_t i = 0; i < 4; i++) {
          TaskFrame task_i;
          task_i.node_1 = node_1_children[i];
          task_i.node_2 = node_2;
          task_i.interval = obb_interval;
          task_stack.push(std::move(task_i));
        }
      } else {
        // Continue on hm_2
        assert(!hm_2.is_bottom_layer(node_2.layer));
        std::array<HeightMapNode, 4> node_2_children;
        // clang-format off
        node_2_children[0] = HeightMapNode{node_2.layer + 1, Pixel(node_2.pixel.x * 2 + 0, node_2.pixel.y * 2 + 0)};
        node_2_children[1] = HeightMapNode{node_2.layer + 1, Pixel(node_2.pixel.x * 2 + 1, node_2.pixel.y * 2 + 0)};
        node_2_children[2] = HeightMapNode{node_2.layer + 1, Pixel(node_2.pixel.x * 2 + 0, node_2.pixel.y * 2 + 1)};
        node_2_children[3] = HeightMapNode{node_2.layer + 1, Pixel(node_2.pixel.x * 2 + 1, node_2.pixel.y * 2 + 1)};
        // clang-format on

        // Add the tasks
        for (std::size_t i = 0; i < 4; i++) {
          TaskFrame task_i;
          task_i.node_1 = node_1;
          task_i.node_2 = node_2_children[i];
          task_i.interval = obb_interval;
          task_stack.push(std::move(task_i));
        }
      }
    }
  }
}

template <typename S>
void TranslationalDisplacementHeightMapSolver<S>::RunHeightMapPair(
    const HeightMapCollisionGeometry<S>* hm1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement,
    const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request_in,
    ContinuousCollisionResult<S>& result_in) {
  // Check input
  if (hm1 == nullptr || hm1->raw_heightmap() == nullptr || hm2 == nullptr ||
      hm2->raw_heightmap() == nullptr) {
    return;
  }

  // Assign data
  request = &request_in;
  result = &result_in;
  if (result->TerminationConditionSatisfied(request->num_max_contacts)) {
    return;
  }

  // Into interface
  runHeightMapPair(hm1, tf1, hm1_displacement, hm2, tf2);
}

}  // namespace detail
}  // namespace fcl