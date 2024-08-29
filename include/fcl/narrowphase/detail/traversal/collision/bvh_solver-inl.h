#pragma once

namespace fcl {
namespace detail {

template <typename BV>
template <typename Shape>
void OrientedNodeBVHSolver<BV>::MeshShapeIntersect(
    const BVHModel<BV>* bvh_1, const Shape& shape_2, const Transform3<S>& tf1,
    const Transform3<S>& tf2, const CollisionRequest<S>& request,
    CollisionResult<S>& result) const {
  if (request.terminationConditionSatisfied(result)) {
    return;
  }

  // Compute the bv for shape_2
  BV shape_bv;
  computeBV(shape_2, tf2, shape_bv);

  // The task stack contains only the bv id
  std::stack<int> bv_id_stack;
  bv_id_stack.push(0);
  ShapePairIntersectSolver<S> shape_solver(gjk_solver);

  // Processing loop
  while (!bv_id_stack.empty()) {
    // Current task
    const int bv_id = bv_id_stack.top();
    bv_id_stack.pop();

    // Check bv for both leaf and internal
    assert(bv_id < bvh_1->getNumBVs());
    const BVNode<BV>& bv_node = bvh_1->getBV(bv_id);
    const BV& node1_bv = bv_node.bv;
    const bool bv_overlap =
        overlap(tf1.linear(), tf1.translation(), shape_bv, node1_bv);
    if (!bv_overlap) {
      continue;
    }

    // Not leaf, continue down
    if (!bv_node.isLeaf()) {
      const int child_1 = bv_node.leftChild();
      const int child_2 = bv_node.rightChild();
      bv_id_stack.push(child_1);
      bv_id_stack.push(child_2);
      continue;
    }

    // Leaf case
    assert(bv_node.isLeaf());
    const int primitive_id = bv_node.primitiveId();
    const Simplex<S> simplex = bvh_1->getSimplex(primitive_id);

    // Make the shape solver functor
    ContactMeta<S> contact;
    contact.o1 = bvh_1;
    contact.o2 = &shape_2;
    contact.b1 = primitive_id;
    contact.b2 = Contact<S>::NONE;
    contact.reverse_normal = true;

    // Go to shape solver
    shape_solver.template ShapeSimplexIntersect<Shape>(
        shape_2, tf2, simplex, tf1, request, contact, result);

    // Check termination
    if (request.terminationConditionSatisfied(result)) {
      return;
    }
  }
}

template <typename BV>
void OrientedNodeBVHSolver<BV>::MeshIntersect(
    const BVHModel<BV>* bvh_1, const BVHModel<BV>* bvh_2,
    const Transform3<S>& tf1, const Transform3<S>& tf2,
    const CollisionRequest<S>& request, CollisionResult<S>& result) const {
  if (request.terminationConditionSatisfied(result)) {
    return;
  }

  // Check initial isOccupied
  if ((!bvh_1->isOccupied()) || (!bvh_2->isOccupied())) {
    return;
  }

  // Compute the relative transform
  Matrix3<S> rotation_2to1;
  Vector3<S> translation_2in1;
  relativeTransform(tf1.linear(), tf1.translation(), tf2.linear(),
                    tf2.translation(), rotation_2to1, translation_2in1);

  // The task stack contains the bv id for both side
  std::stack<std::pair<int, int>> bv_id_pair_stack;
  bv_id_pair_stack.emplace(0, 0);
  ShapePairIntersectSolver<S> shape_solver(gjk_solver);

  // Processing loop
  while (!bv_id_pair_stack.empty()) {
    // Gather the id
    const auto& task_pair = bv_id_pair_stack.top();
    const int bv_id_1 = task_pair.first;
    const int bv_id_2 = task_pair.second;
    bv_id_pair_stack.pop();

    // Gather the node
    assert(bv_id_1 < bvh_1->getNumBVs());
    assert(bv_id_2 < bvh_2->getNumBVs());
    const BVNode<BV>& bv_node_1 = bvh_1->getBV(bv_id_1);
    const BVNode<BV>& bv_node_2 = bvh_2->getBV(bv_id_2);

    // Check bv
    const BV& bv_1 = bv_node_1.bv;
    const BV& bv_2 = bv_node_2.bv;
    const bool bv_overlap =
        overlap(rotation_2to1, translation_2in1, bv_1, bv_2);
    if (!bv_overlap) {
      continue;
    }

    // Need to go down
    const bool is_leaf_1 = bv_node_1.isLeaf();
    const bool is_leaf_2 = bv_node_2.isLeaf();
    if (is_leaf_1 && is_leaf_2) {
      const int primitive_id1 = bv_node_1.primitiveId();
      const int primitive_id2 = bv_node_2.primitiveId();
      const Simplex<S>& simplex1 = bvh_1->getSimplex(primitive_id1);
      const Simplex<S>& simplex2 = bvh_2->getSimplex(primitive_id2);

      // Make the shape solver
      ContactMeta<S> contact;
      contact.o1 = bvh_1;
      contact.o2 = bvh_2;
      contact.b1 = primitive_id1;
      contact.b2 = primitive_id2;

      // Go to gjk solver
      shape_solver.SimplexIntersect(simplex1, tf1, simplex2, tf2, rotation_2to1,
                                    translation_2in1, request, contact,
                                    result);

      // Check termination
      if (request.terminationConditionSatisfied(result)) {
        return;
      }
    } else {
      // At least one is not leaf, continue
      const bool continue_on_1 =
          (is_leaf_2) || ((!is_leaf_1) && bv_1.size() > bv_2.size());
      if (continue_on_1) {
        const int child_1 = bv_node_1.leftChild();
        const int child_2 = bv_node_1.rightChild();
        bv_id_pair_stack.push(std::make_pair(child_1, bv_id_2));
        bv_id_pair_stack.push(std::make_pair(child_2, bv_id_2));
      } else {
        const int child_1 = bv_node_2.leftChild();
        const int child_2 = bv_node_2.rightChild();
        bv_id_pair_stack.push(std::make_pair(bv_id_1, child_1));
        bv_id_pair_stack.push(std::make_pair(bv_id_1, child_2));
      }
    }
  }
}

}  // namespace detail
}  // namespace fcl