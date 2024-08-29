//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
std::int64_t CollisionSolverOctree2<S>::encodeOctree2Node(
    std::uint32_t octree_node_vector_index, bool is_leaf_node,
    std::uint8_t inner_child_idx) {
  std::int64_t encoded = octree_node_vector_index;
  encoded += (static_cast<std::int64_t>(inner_child_idx) << 32);
  const std::uint8_t is_leaf_node_int = is_leaf_node ? 1 : 0;
  encoded += (static_cast<std::int64_t>(is_leaf_node_int) << 48);
  return encoded;
}

template <typename S>
template <typename Shape>
void CollisionSolverOctree2<S>::boxToShapeProcessLeafPair(
    const Octree2CollisionGeometry<S>& octree, const Transform3<S>& tf_octree,
    const Shape& shape, const Transform3<S>& tf_shape,
    std::int64_t encoded_octree_node_idx, const AABB<S>& voxel_aabb,
    OctreeLeafComputeCache& cache) const {
  Box<S>& box = cache.box1;
  Transform3<S>& box_tf = cache.shape1_tf;
  constructBox(voxel_aabb, tf_octree, box, box_tf);

  // Run solver
  auto& contact = cache.contact_meta;
  contact.reset();
  contact.o1 = &octree;
  contact.o2 = &shape;
  contact.b1 = encoded_octree_node_idx;
  contact.b2 = Contact<S>::NONE;
  contact.o1_bv = voxel_aabb;
  cache.shape_solver.template ShapeIntersect<Box<S>, Shape>(
      box, box_tf, shape, tf_shape, *request, contact, *result);
}

template <typename S>
void CollisionSolverOctree2<S>::boxToSimplexProcessLeafPair(
    const Octree2CollisionGeometry<S>& octree, const Transform3<S>& tf_octree,
    const AABB<S>& voxel_aabb, std::int64_t encoded_octree_node_idx,
    const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
    const Simplex<S>& simplex, fcl::intptr_t index_2,
    OctreeLeafComputeCache& cache) const {
  Box<S>& box = cache.box1;
  Transform3<S>& box_tf = cache.shape1_tf;
  constructBox(voxel_aabb, tf_octree, box, box_tf);

  // Run solver
  auto& contact = cache.contact_meta;
  contact.reset();
  contact.o1 = &octree;
  contact.o2 = bvh;
  contact.b1 = encoded_octree_node_idx;
  contact.b2 = index_2;
  contact.o1_bv = voxel_aabb;
  cache.shape_solver.template ShapeSimplexIntersect<Box<S>>(
      box, box_tf, simplex, tf_bvh, *request, contact, *result);
}

template <typename S>
bool sphereApproximationVoxelOverlap(const AABB<S>& bv_1, const AABB<S>& bv_2,
                                     const Transform3<S>& tf1,
                                     const Transform3<S>& tf2) {
  const S radius_1 = S(0.886) * (bv_1.max_[0] - bv_1.min_[0]);
  const S radius_2 = S(0.886) * (bv_2.max_[0] - bv_2.min_[0]);
  const Vector3<S> center_diff = tf1 * bv_1.center() - tf2 * bv_2.center();
  const S distance_square = center_diff.x() * center_diff.x() +
                            center_diff.y() * center_diff.y() +
                            center_diff.z() * center_diff.z();
  bool sphere_overlap_out =
      distance_square < (radius_1 + radius_2) * (radius_1 + radius_2);
  return sphere_overlap_out;
}

template <typename S>
bool CollisionSolverOctree2<S>::octreePairTwoLeafNode(
    const Octree2CollisionGeometry<S>& tree1,
    const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
    const Transform3<S>& tf_tree2,
    const octree2::OctreeTraverseStackElement<S>& tree1_elem,
    const octree2::OctreeTraverseStackElement<S>& tree2_elem,
    const octree2::OctreeLeafNode& tree1_node,
    const octree2::OctreeLeafNode& tree2_node,
    const FixedRotationBoxDisjoint<S>& disjoint,
    OctreeLeafComputeCache& cache) const {
  // Make the shape solver
  Transform3<S>& shape1_tf = cache.shape1_tf;
  Transform3<S>& shape2_tf = cache.shape2_tf;
  Box<S>& box1 = cache.box1;
  Box<S>& box2 = cache.box2;
  const bool request_penetration = request->isPenetrationEnabled();

  // Meta fn
  auto& contact = cache.contact_meta;
  contact.reset();
  contact.o1 = &tree1;
  contact.o2 = &tree2;
  AABB<S>& o1_bv = contact.o1_bv;
  AABB<S>& o2_bv = contact.o2_bv;

  if (tree1_node.is_fully_occupied() && tree2_node.is_fully_occupied()) {
    o1_bv = tree1_elem.bv;
    o2_bv = tree2_elem.bv;
    contact.b1 = encodeOctree2Node(tree1_elem.node_vector_index, true);
    contact.b2 = encodeOctree2Node(tree2_elem.node_vector_index, true);

    // Check whether penetration is required
    if (!request_penetration) {
      if (result->numContacts() < request->maxNumContacts() &&
          (!disjoint.isDisjoint(o1_bv, o2_bv, true))) {
        Contact<S> this_contact;
        contact.writeToContact(this_contact);
        result->addContact(this_contact);
      }
    } else {
      // Make box and run solver
      constructBox(o1_bv, tf_tree1, box1, shape1_tf);
      constructBox(o2_bv, tf_tree2, box2, shape2_tf);
      cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
          box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
    }

    // Done with this case
    return request->terminationConditionSatisfied(*result);
  }

  if (tree1_node.is_fully_occupied() && (!tree2_node.is_fully_occupied())) {
    // Assign the bv for o1
    o1_bv = tree1_elem.bv;
    contact.b1 = encodeOctree2Node(tree1_elem.node_vector_index, true);

    // Make the shape for o1
    if (request_penetration) constructBox(o1_bv, tf_tree1, box1, shape1_tf);

    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      if (!tree2_node.child_occupied.test_i(child_i)) continue;
      octree2::computeChildAABB(tree2_elem.bv, child_i, o2_bv);
      contact.b2 =
          encodeOctree2Node(tree2_elem.node_vector_index, true, child_i);

      if (request_penetration) {
        constructBox(o2_bv, tf_tree2, box2, shape2_tf);
        cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
            box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
      } else {
        if (result->numContacts() < request->maxNumContacts() &&
            (!disjoint.isDisjoint(o1_bv, o2_bv, true))) {
          Contact<S> this_contact;
          contact.writeToContact(this_contact);
          result->addContact(this_contact);
        }
      }

      // Done with this case
      if (request->terminationConditionSatisfied(*result)) return true;
    }

    // Not completely finished
    return false;
  }

  if (tree2_node.is_fully_occupied() && (!tree1_node.is_fully_occupied())) {
    // Assign the bv for o2
    o2_bv = tree2_elem.bv;
    contact.b2 = encodeOctree2Node(tree2_elem.node_vector_index, true);
    if (request_penetration) constructBox(o2_bv, tf_tree2, box2, shape2_tf);

    // Iterate the children
    for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
      if (!tree1_node.child_occupied.test_i(child_i)) continue;
      octree2::computeChildAABB(tree1_elem.bv, child_i, o1_bv);
      contact.b1 =
          encodeOctree2Node(tree1_elem.node_vector_index, true, child_i);

      if (request_penetration) {
        constructBox(o1_bv, tf_tree1, box1, shape1_tf);
        cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
            box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
      } else {
        if (result->numContacts() < request->maxNumContacts() &&
            (!disjoint.isDisjoint(o1_bv, o2_bv, true))) {
          Contact<S> this_contact;
          contact.writeToContact(this_contact);
          result->addContact(this_contact);
        }
      }

      // Done with this case
      if (request->terminationConditionSatisfied(*result)) return true;
    }

    // Not completely finished
    return false;
  }

  // Both are not fully satisfied
  assert(!tree1_node.is_fully_occupied());
  assert(!tree2_node.is_fully_occupied());
  if (disjoint.isDisjoint(tree1_elem.bv, tree2_elem.bv, false)) return false;

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
        const bool overlap = sphereApproximationVoxelOverlap(
            o1_bv, tree2_elem.bv, tf_tree1, tf_tree2);
        if (!overlap) tree1_child2check.clear_i(i);
      }

      // Prune tree2 child
      for (std::uint8_t i = 0; i < 8; i++) {
        if (!tree2_child2check.test_i(i)) continue;
        octree2::computeChildAABB(tree2_elem.bv, i, o2_bv);
        const bool overlap = sphereApproximationVoxelOverlap(
            tree1_elem.bv, o2_bv, tf_tree1, tf_tree2);
        if (!overlap) tree2_child2check.clear_i(i);
      }
    }
  }

  // Loop over child
  for (std::uint8_t child_1 = 0; child_1 < 8; child_1++) {
    if (!tree1_child2check.test_i(child_1)) continue;
    octree2::computeChildAABB(tree1_elem.bv, child_1, o1_bv);
    contact.b1 = encodeOctree2Node(tree1_elem.node_vector_index, true, child_1);
    if (request_penetration) constructBox(o1_bv, tf_tree1, box1, shape1_tf);

    // Into child 2
    for (std::uint8_t child_2 = 0; child_2 < 8; child_2++) {
      if (!tree2_child2check.test_i(child_2)) continue;
      octree2::computeChildAABB(tree2_elem.bv, child_2, o2_bv);
      contact.b2 =
          encodeOctree2Node(tree2_elem.node_vector_index, true, child_2);

      if (request_penetration) {
        constructBox(o2_bv, tf_tree2, box2, shape2_tf);
        cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
            box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
      } else {
        if (result->numContacts() < request->maxNumContacts() &&
            (!disjoint.isDisjoint(o1_bv, o2_bv, true))) {
          Contact<S> this_contact;
          contact.writeToContact(this_contact);
          result->addContact(this_contact);
        }
      }

      // Done after update
      if (request->terminationConditionSatisfied(*result)) return true;
    }
  }

  // Done
  return false;
}

template <typename S>
bool CollisionSolverOctree2<S>::octreePairInnerNodeWithLeafNode(
    const Octree2CollisionGeometry<S>& tree1,
    const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
    const Transform3<S>& tf_tree2,
    const octree2::OctreeTraverseStackElement<S>& tree1_elem,
    const octree2::OctreeTraverseStackElement<S>& tree2_elem,
    const octree2::OctreeLeafNode& tree2_node, bool reverse_tree12,
    const FixedRotationBoxDisjoint<S>& disjoint,
    OctreeLeafComputeCache& cache) const {
  Transform3<S>& shape1_tf = cache.shape1_tf;
  Transform3<S>& shape2_tf = cache.shape2_tf;
  Box<S>& box1 = cache.box1;
  Box<S>& box2 = cache.box2;

  // Meta fn
  AABB<S>& leaf2_bv = cache.aabb2;
  auto& contact = cache.contact_meta;
  contact.reset();
  contact.o1 = &tree1;
  contact.o2 = &tree2;
  contact.b1 = encodeOctree2Node(tree1_elem.node_vector_index, false);
  contact.o1_bv = tree1_elem.bv;
  contact.reverse_o1_and_o2 = reverse_tree12;
  contact.reverse_normal = reverse_tree12;

  // Make box or obb, by request
  const bool request_penetration = request->isPenetrationEnabled();
  if (request_penetration)
    constructBox(tree1_elem.bv, tf_tree1, box1, shape1_tf);

  // If this is a full node
  if (tree2_node.is_fully_occupied()) {
    // Assign the bv
    leaf2_bv = tree2_elem.bv;
    contact.o2_bv = leaf2_bv;
    contact.b2 = encodeOctree2Node(tree2_elem.node_vector_index, true);

    if (request_penetration) {
      constructBox(leaf2_bv, tf_tree2, box2, shape2_tf);
      cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
          box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
    } else {
      const bool is_disjoint =
          reverse_tree12 ? disjoint.isDisjoint(leaf2_bv, tree1_elem.bv, true)
                         : disjoint.isDisjoint(tree1_elem.bv, leaf2_bv, true);
      if (result->numContacts() < request->maxNumContacts() && (!is_disjoint)) {
        Contact<S> this_contact;
        contact.writeToContact(this_contact);
        result->addContact(this_contact);
      }
    }

    // Done with this case
    return request->terminationConditionSatisfied(*result);
  }

  // Iterate over nodes in tree2
  assert(!tree2_node.is_fully_occupied());
  for (std::uint8_t child_i = 0; child_i < 8; child_i++) {
    // Check occupied
    if (!tree2_node.child_occupied.test_i(child_i)) continue;

    // Always compute the child AABB
    octree2::computeChildAABB(tree2_elem.bv, child_i, leaf2_bv);
    contact.o2_bv = leaf2_bv;
    contact.b2 = encodeOctree2Node(tree2_elem.node_vector_index, true, child_i);

    if (request_penetration) {
      constructBox(leaf2_bv, tf_tree2, box2, shape2_tf);
      cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
          box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
    } else {
      const bool is_disjoint =
          reverse_tree12 ? disjoint.isDisjoint(leaf2_bv, tree1_elem.bv, true)
                         : disjoint.isDisjoint(tree1_elem.bv, leaf2_bv, true);
      if (result->numContacts() < request->maxNumContacts() && (!is_disjoint)) {
        Contact<S> this_contact;
        contact.writeToContact(this_contact);
        result->addContact(this_contact);
      }
    }

    // Done after update
    if (request->terminationConditionSatisfied(*result)) return true;
  }

  // Not done
  return false;
}

template <typename S>
void CollisionSolverOctree2<S>::octreePairInnerNodePairAsLeaf(
    const Octree2CollisionGeometry<S>& tree1,
    const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
    const Transform3<S>& tf_tree2,
    const octree2::OctreeTraverseStackElement<S>& tree1_elem,
    const octree2::OctreeTraverseStackElement<S>& tree2_elem,
    const FixedRotationBoxDisjoint<S>& disjoint,
    OctreeLeafComputeCache& cache) const {
  // If sufficient count, this function does nothing
  if (result->numContacts() >= request->maxNumContacts()) return;

  // Contact meta
  auto& contact = cache.contact_meta;
  contact.reset();
  contact.o1 = &tree1;
  contact.o2 = &tree2;
  contact.b1 = encodeOctree2Node(tree1_elem.node_vector_index, false);
  contact.b2 = encodeOctree2Node(tree2_elem.node_vector_index, false);
  contact.o1_bv = tree1_elem.bv;
  contact.o2_bv = tree2_elem.bv;

  if (request->isPenetrationEnabled()) {
    Transform3<S>& shape1_tf = cache.shape1_tf;
    Transform3<S>& shape2_tf = cache.shape2_tf;
    Box<S>& box1 = cache.box1;
    Box<S>& box2 = cache.box2;
    constructBox(tree1_elem.bv, tf_tree1, box1, shape1_tf);
    constructBox(tree2_elem.bv, tf_tree2, box2, shape2_tf);
    cache.shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
        box1, shape1_tf, box2, shape2_tf, *request, contact, *result);
  } else {
    if (!disjoint.isDisjoint(tree1_elem.bv, tree2_elem.bv, true)) {
      Contact<S> this_contact;
      contact.writeToContact(this_contact);
      result->addContact(this_contact);
    }
  }
}

}  // namespace detail
}  // namespace fcl
