#pragma once

namespace fcl {
namespace detail {

struct HeightMapNode {
  int layer{0};
  heightmap::Pixel pixel;
  HeightMapNode() = default;
  HeightMapNode(int layer_in, heightmap::Pixel pixel_in)
      : layer(layer_in), pixel(pixel_in) {}
};

template <typename S>
bool heightMapNodeToAABB(const heightmap::LayeredHeightMap<S>& height_map,
                         const HeightMapNode& node, AABB<S>& aabb) {
  const heightmap::FlatHeightMap<S>& hm_layer = height_map.layers()[node.layer];
  return hm_layer.pixelToBox(node.pixel, aabb);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::heightMapShapeIntersectImpl(
    const HeightMapCollisionGeometry<S>& heightmap_geometry, const Shape& shape,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
    ShapeHeightMapCollisionType option_code) const {
  // Construct bv for shape
  const heightmap::FlatHeightMap<S>& heightmap =
      heightmap_geometry.raw_heightmap()->bottom();
  Transform3<S> tf_shape_to_map = tf_hm.inverse() * tf_shape;
  AABB<S> shape_aabb_in_hm;
  computeBV(shape, tf_shape_to_map, shape_aabb_in_hm);

  // Check z
  if (shape_aabb_in_hm.max_.z() < 0) return;
  if (shape_aabb_in_hm.min_.z() > heightmap.height_upper_bound_meter()) return;

  // Check roi and visit
  heightmap::PixelSpaceROI roi;
  bool has_overlap_2d = heightmap.aabbToPixelROI(shape_aabb_in_hm, roi);
  if (!has_overlap_2d) return;

  // Determine which one to use
  bool use_flat_impl = (option_code == ShapeHeightMapCollisionType::Flat);
  if ((!use_flat_impl) &&
      (option_code != ShapeHeightMapCollisionType::Traversal)) {
    // Determine automatically by roi size
    std::size_t total_size =
        heightmap.full_shape_x() * heightmap.full_shape_y();
    uint16_t roi_delta_x = roi.bottom_right.x >= roi.top_left.x
                               ? (roi.bottom_right.x - roi.top_left.x + 1)
                               : 0;
    uint16_t roi_delta_y = roi.bottom_right.y >= roi.top_left.y
                               ? (roi.bottom_right.y - roi.top_left.y + 1)
                               : 0;
    std::size_t roi_size = roi_delta_x * roi_delta_y;
    assert(roi_size > 0);

    constexpr std::size_t kSmallROI_Ratio = 8;
    if (roi_size * kSmallROI_Ratio < total_size) {
      use_flat_impl = true;
    }
  }

  // Forward to impl
  if (use_flat_impl) {
    flatHeightMapShapeIntersectImpl(heightmap_geometry, shape, shape_aabb_in_hm,
                                    roi, tf_hm, tf_shape);
  } else {
    traversalHeightMapShapeIntersectImpl(
        heightmap_geometry, shape, shape_aabb_in_hm, roi, tf_hm, tf_shape);
  }
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::flatHeightMapShapeIntersectImpl(
    const HeightMapCollisionGeometry<S>& heightmap_geometry, const Shape& shape,
    const AABB<S>& shape_aabb_in_hm, const heightmap::PixelSpaceROI& bottom_roi,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape) const {
  // Visit functor
  const S half_resolution_x =
      S(0.5) * heightmap_geometry.raw_heightmap()->bottom().resolution_x();
  const S half_resolution_y =
      S(0.5) * heightmap_geometry.raw_heightmap()->bottom().resolution_y();
  auto check_collision = [&](const heightmap::Pixel& pixel,
                             const heightmap::Point2D<S>& box_bottom_center,
                             uint16_t height_in_mm) -> bool {
    // Basic check
    (void)(box_bottom_center);
    if (height_in_mm == 0) return false;
    const S height_in_meter = static_cast<S>(height_in_mm) * S(0.001);
    if (shape_aabb_in_hm.min_.z() > height_in_meter) {
      return false;
    }

    // Compute AABB
    AABB<S> pixel_aabb;
    pixel_aabb.min_.x() = box_bottom_center.x() - half_resolution_x;
    pixel_aabb.max_.x() = box_bottom_center.x() + half_resolution_x;
    pixel_aabb.min_.y() = box_bottom_center.y() - half_resolution_y;
    pixel_aabb.max_.y() = box_bottom_center.y() + half_resolution_y;
    pixel_aabb.min_.z() = S(0.0);
    pixel_aabb.max_.z() = height_in_meter;

    // Process leaves
    boxToShapeProcessLeafPair(heightmap_geometry, tf_hm, shape, tf_shape, pixel,
                              pixel_aabb);
    return request->terminationConditionSatisfied(*result);
  };

  // Visit within roi
  heightmap_geometry.raw_heightmap()->bottom().visitHeightMap(check_collision,
                                                              &bottom_roi);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::traversalHeightMapShapeIntersectImpl(
    const HeightMapCollisionGeometry<S>& heightmap_geometry, const Shape& shape,
    const AABB<S>& shape_aabb_in_hm, const heightmap::PixelSpaceROI& bottom_roi,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape) const {
  // The task stack contains the element for checking
  (void)(bottom_roi);
  // The task stack contains the element for checking
  using heightmap::Pixel;
  const LayeredHeightMap& heightmap = *heightmap_geometry.raw_heightmap();
  std::stack<HeightMapNode> task_stack;

  // Init the task stack
  {
    const FlatHeightMap& hm_top = heightmap.top();
    for (uint16_t y = 0; y < hm_top.full_shape_y(); y++) {
      for (uint16_t x = 0; x < hm_top.full_shape_x(); x++) {
        task_stack.push(HeightMapNode(0, Pixel(x, y)));
      }
    }
  }

  // Compute obb
  FixedRotationBoxDisjoint<S> obb_disjoint;
  AABB<S> shape_obb_local_AABB;
  {
    OBB<S> shape_obb_world;
    computeBV(shape, tf_shape, shape_obb_world);

    // Convert OBB to local AABB and a tf_AABB on that AABB
    shape_obb_local_AABB.max_ = shape_obb_world.extent;
    shape_obb_local_AABB.min_ = -shape_obb_world.extent;

    // Make tf_AABB frame
    Transform3<S> tf_shape_AABB;
    tf_shape_AABB.setIdentity();
    tf_shape_AABB.linear().matrix() = shape_obb_world.axis;
    tf_shape_AABB.translation() = shape_obb_world.To;
    obb_disjoint.initialize(tf_hm, tf_shape_AABB);
  }

  // The processing loop
  while (!task_stack.empty()) {
    // Pop the stack
    const HeightMapNode hm_node = task_stack.top();
    task_stack.pop();
    AABB<S> aabb_local;
    if (!heightMapNodeToAABB<S>(heightmap, hm_node, aabb_local)) continue;
    if (!shape_aabb_in_hm.overlap(aabb_local)) continue;

    // Check OBB
    if (obb_disjoint.isDisjoint(aabb_local, shape_obb_local_AABB, false))
      continue;

    if (heightmap.is_bottom_layer(hm_node.layer)) {
      // Both are leaf nodes
      boxToShapeProcessLeafPair(heightmap_geometry, tf_hm, shape, tf_shape,
                                hm_node.pixel, aabb_local);
      if (request->terminationConditionSatisfied(*result)) return;
    } else {  // at least one is not leaf
      // clang-format off
      task_stack.push(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 0, hm_node.pixel.y * 2 + 0)});
      task_stack.push(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 1, hm_node.pixel.y * 2 + 0)});
      task_stack.push(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 0, hm_node.pixel.y * 2 + 1)});
      task_stack.push(HeightMapNode{hm_node.layer + 1, Pixel(hm_node.pixel.x * 2 + 1, hm_node.pixel.y * 2 + 1)});
      // clang-format on
    }
  }
}

template <typename S>
void HeightMapCollisionSolver<S>::heightMapPairIntersect(
    const HeightMapCollisionGeometry<S>& hm_1_geometry,
    const HeightMapCollisionGeometry<S>& hm_2_geometry,
    const Transform3<S>& tf1, const Transform3<S>& tf2) const {
  // The task stack contains the element for checking
  using heightmap::Pixel;
  const LayeredHeightMap& hm_1 = *hm_1_geometry.raw_heightmap();
  const LayeredHeightMap& hm_2 = *hm_2_geometry.raw_heightmap();
  struct TaskFrame {
    HeightMapNode node_1;
    HeightMapNode node_2;
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
            task_stack.push(std::move(this_task));
          }
        }
      }
    }
  }

  // Init for disjoint
  FixedRotationBoxDisjoint<S> obb_disjoint;
  obb_disjoint.initialize(tf1, tf2);

  // The processing loop
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskFrame this_task = task_stack.top();
    task_stack.pop();
    const HeightMapNode node_1 = this_task.node_1;
    const HeightMapNode node_2 = this_task.node_2;
    AABB<S> aabb_1, aabb_2;
    const bool aabb_1_valid = heightMapNodeToAABB<S>(hm_1, node_1, aabb_1);
    const bool aabb_2_valid = heightMapNodeToAABB<S>(hm_2, node_2, aabb_2);
    if ((!aabb_1_valid) || (!aabb_2_valid)) continue;

    if (hm_1.is_bottom_layer(node_1.layer) &&
        hm_2.is_bottom_layer(node_2.layer)) {
      // Both are leaf nodes
      boxToBoxProcessLeafPair(
          &hm_1_geometry, tf1, aabb_1, heightmap::encodePixel(node_1.pixel),
          &hm_2_geometry, tf2, aabb_2, heightmap::encodePixel(node_2.pixel),
          obb_disjoint);
      if (request->terminationConditionSatisfied(*result)) return;
    } else {  // at least one is not leaf
      // Check OBB bv
      if (obb_disjoint.isDisjoint(aabb_1, aabb_2, false)) continue;

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
          task_stack.push(std::move(task_i));
        }
      }
    }

    // The stack loop
  }
}

template <typename S>
template <typename BV>
void HeightMapCollisionSolver<S>::heightMapBVHIntersect(
    const HeightMapCollisionGeometry<S>* hm_geometry, const BVHModel<BV>* bvh,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_bvh) const {
  assert(hm_geometry != nullptr);
  assert(bvh != nullptr);

  // We need a stack following the recursive procedure
  struct TaskFrame {
    HeightMapNode hm_node;
    int bvh_node;
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
        task_stack.push(std::move(this_task));
      }
    }
  }

  // Stack loop
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskFrame this_task = task_stack.top();
    task_stack.pop();
    const HeightMapNode hm_node_1 = this_task.hm_node;
    const int bvh_node_2 = this_task.bvh_node;

    // Compute the bv
    AABB<S> hm_aabb_1;
    const bool hm_aabb_1_valid =
        heightMapNodeToAABB<S>(hm_1, hm_node_1, hm_aabb_1);
    if (!hm_aabb_1_valid) continue;
    const BV node_2_bv = bvh->getBV(bvh_node_2).bv;

    // Compute obb
    OBB<S> obb_1_hm, obb_2_bvh;
    convertBV(hm_aabb_1, tf_hm, obb_1_hm);
    convertBV(node_2_bv, tf_bvh, obb_2_bvh);

    // Always check bvh
    if (!obb_1_hm.overlap(obb_2_bvh)) {
      continue;
    }

    // The root case
    if (bvh->getBV(bvh_node_2).isLeaf() &&
        hm_1.is_bottom_layer(hm_node_1.layer)) {
      // Obb overlap checked. Just obtain the triangle is supported
      const int primitive_id = bvh->getBV(bvh_node_2).primitiveId();
      const Simplex<S> simplex = bvh->getSimplex(primitive_id);
      boxToTriangleProcessLeafPair(hm_geometry, tf_hm, hm_aabb_1,
                                   heightmap::encodePixel(hm_node_1.pixel), bvh,
                                   tf_bvh, simplex, primitive_id);

      // Finished
      if (request->terminationConditionSatisfied(*result)) return;
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
          task_stack.push(std::move(task_i));
        }
      } else {
        // Continue on bvh_2, push left child to stack
        TaskFrame new_task_left;
        new_task_left.hm_node = hm_node_1;
        new_task_left.bvh_node = bvh->getBV(bvh_node_2).leftChild();
        task_stack.push(std::move(new_task_left));

        // Push right child to stack
        TaskFrame new_task_right;
        new_task_right.hm_node = hm_node_1;
        new_task_right.bvh_node = bvh->getBV(bvh_node_2).rightChild();
        task_stack.push(std::move(new_task_right));
      }
    }

    // The stack loop
  }
}

template <typename S>
void HeightMapCollisionSolver<S>::heightMapOctreeIntersect(
    const HeightMapCollisionGeometry<S>* hm_geometry,
    const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_hm,
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
        task_stack.push(std::move(this_task));
      }
    }
  }

  // Make OBB disjoint
  FixedRotationBoxDisjoint<S> obb_disjoint;
  obb_disjoint.initialize(tf_hm, tf_octree);

  // Task processing loop
  AABB<S> aabb_1_hm, child_aabb_tree_2;
  while (!task_stack.empty()) {
    // Pop the stack
    const TaskStackElement this_task = task_stack.top();
    task_stack.pop();

    // Obtain the task data
    const auto hm_node_1 = this_task.hm_node_1;
    const auto octree_node_2 = this_task.octree_node_2;

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
      // This is inner node aas leaf
      if (!octree_node_2.is_leaf_node) {
        assert(inner_nodes_full[octree_node_2.node_vector_index]);
        boxToBoxProcessLeafPair(hm_geometry, tf_hm, aabb_1_hm,
                                heightmap::encodePixel(hm_node_1.pixel), octree,
                                tf_octree, octree_node_2.bv,
                                octree_node_2.node_vector_index, obb_disjoint);

        if (request->terminationConditionSatisfied(*result)) return;
        continue;  // finish this case
      }

      // This is real leaf
      assert(octree_node_2.is_leaf_node);
      const auto& leaf_2 = leaf_nodes[octree_node_2.node_vector_index];
      if (leaf_2.is_fully_occupied()) {
        boxToBoxProcessLeafPair(hm_geometry, tf_hm, aabb_1_hm,
                                heightmap::encodePixel(hm_node_1.pixel), octree,
                                tf_octree, octree_node_2.bv,
                                octree_node_2.node_vector_index, obb_disjoint);

        if (request->terminationConditionSatisfied(*result)) return;
      } else {
        for (std::uint8_t i = 0; i < 8; i++) {
          if (!leaf_2.child_occupied.test_i(i)) continue;
          octree2::computeChildAABB(octree_node_2.bv, i, child_aabb_tree_2);
          boxToBoxProcessLeafPair(
              hm_geometry, tf_hm, aabb_1_hm,
              heightmap::encodePixel(hm_node_1.pixel), octree, tf_octree,
              child_aabb_tree_2, octree_node_2.node_vector_index, obb_disjoint);
          if (request->terminationConditionSatisfied(*result)) return;
        }
      }

      // Finish the leaf case
      continue;
    }

    // Non-leaf case
    const bool is_obb_disjoint =
        obb_disjoint.isDisjoint(aabb_1_hm, octree_node_2.bv, true);
    if (is_obb_disjoint) continue;

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
      task_stack.push(std::move(new_stack_element));
    }
  }
}

}  // namespace detail
}  // namespace fcl