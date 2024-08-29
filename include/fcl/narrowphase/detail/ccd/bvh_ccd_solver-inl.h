#pragma once

#include "fcl/narrowphase/detail/ccd/ccd_solver_utility.h"

namespace fcl {
namespace detail {

/// For shape-mesh
template <typename S, typename Shape, typename BV>
struct TranslationalDisplacementShapeBVHSolverImpl {
  // Does not support general bv
  static void Run(const Shape*, const Transform3<S>&,
                  const TranslationalDisplacement<S>&, const BVHModel<BV>*,
                  const Transform3<S>&, const ContinuousCollisionRequest<S>&,
                  ContinuousCollisionResult<S>&) {
    static_assert(std::is_same<typename Shape::S, S>::value,
                  "scalar type must match");
    std::cerr << "Not supported" << std::endl;
  }
};

template <typename S, typename Shape>
struct TranslationalDisplacementShapeBVHSolverImpl<S, Shape, AABB<S>> {
  // Does not support general bv
  static void Run(const Shape* s1, const Transform3<S>& tf1,
                  const TranslationalDisplacement<S>& s1_displacement,
                  const BVHModel<AABB<S>>* bvh_2, const Transform3<S>& tf2,
                  const ContinuousCollisionRequest<S>& request,
                  ContinuousCollisionResult<S>& result) {
    if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
      return;
    }

    // Init the shape AABB and disjoint
    FixedOrientationBoxPairTranslationalCCD<S> disjoint;
    AABB<S> shape_local_AABB;
    initializeShapeFixedOrientationBoxTranslationalCCD<Shape>(
        *s1, tf1, s1_displacement, tf2, disjoint, shape_local_AABB);

    // The task stack contains bv id and checked interval
    std::stack<std::pair<int, Interval<S>>> bv_id_interval_stack;
    {
      std::pair<int, Interval<S>> init_elem;
      init_elem.first = 0;
      init_elem.second.lower_bound = 0.0;
      init_elem.second.upper_bound = 1.0;
      bv_id_interval_stack.push(std::move(init_elem));
    }

    // Processing loop
    while (!bv_id_interval_stack.empty()) {
      // Current task
      int bv_id = -1;
      Interval<S> bv_parent_interval;
      {
        const auto& stack_top = bv_id_interval_stack.top();
        bv_id = stack_top.first;
        bv_parent_interval = stack_top.second;
        bv_id_interval_stack.pop();
      }

      // Check bv for both leaf and internal
      assert(bv_id < bvh_2->getNumBVs());
      const BVNode<AABB<S>>& bv2_node = bvh_2->getBV(bv_id);

      // Get node2 bv and transform with tf
      const AABB<S>& node2_bv = bv2_node.bv;
      Interval<S> obb_interval;
      const bool is_disjoint =
          disjoint.IsDisjoint(shape_local_AABB, node2_bv, bv_parent_interval,
                              obb_interval, request.zero_movement_tolerance);
      if (is_disjoint) {
        continue;
      }

      // Not leaf, continue down
      if (!bv2_node.isLeaf()) {
        const int child_1 = bv2_node.leftChild();
        const int child_2 = bv2_node.rightChild();
        bv_id_interval_stack.push({child_1, obb_interval});
        bv_id_interval_stack.push({child_2, obb_interval});
        continue;
      }

      // Leaf case
      assert(bv2_node.isLeaf());
      const int primitive_id = bv2_node.primitiveId();
      const Simplex<S>& simplex = bvh_2->getSimplex(primitive_id);

      // Make the shape solver functor
      ContinuousContactMeta<S> contact;
      contact.o1 = s1;
      contact.o2 = bvh_2;
      contact.b1 = ContinuousContactMeta<S>::kNone;
      contact.b2 = primitive_id;
      contact.external_box_toc = obb_interval;

      // Invoke caller
      ShapePairTranslationalCollisionSolver<S>::template RunShapeSimplex<Shape>(
          *s1, tf1, s1_displacement, simplex, tf2, contact, request, result);

      // Check termination
      if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
        return;
      }
    }
  }
};

template <typename S, typename Shape>
struct TranslationalDisplacementShapeBVHSolverImpl<S, Shape, OBB<S>> {
  // Does not support general bv
  static void Run(const Shape* s1, const Transform3<S>& tf1,
                  const TranslationalDisplacement<S>& s1_displacement,
                  const BVHModel<OBB<S>>* bvh_2, const Transform3<S>& tf2,
                  const ContinuousCollisionRequest<S>& request,
                  ContinuousCollisionResult<S>& result) {
    RunSweptBV(s1, tf1, s1_displacement, bvh_2, tf2, request, result);
    // RunExplicitBV(s1, tf1, s1_displacement, bvh_2, tf2, request, result);
  }

  static void RunSweptBV(const Shape* s1, const Transform3<S>& tf1,
                         const TranslationalDisplacement<S>& s1_displacement,
                         const BVHModel<OBB<S>>* bvh_2,
                         const Transform3<S>& tf2,
                         const ContinuousCollisionRequest<S>& request,
                         ContinuousCollisionResult<S>& result) {
    if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
      return;
    }

    // Compute bv and displacement in bv
    OBB<S> shape_bv;
    computeBV(*s1, tf1, shape_bv);

    // Compute the displacement in obb
    TranslationalDisplacement<S> shape_bv_displacement;
    {
      shape_bv_displacement.unit_axis_in_shape1 =
          shape_bv.axis.transpose() *
          (tf1.linear() * s1_displacement.unit_axis_in_shape1);
      shape_bv_displacement.scalar_displacement =
          s1_displacement.scalar_displacement;
    }

    // The task stack contains bv id and checked interval
    std::stack<std::pair<int, Interval<S>>> bv_id_interval_stack;
    {
      std::pair<int, Interval<S>> init_elem;
      init_elem.first = 0;
      init_elem.second.lower_bound = 0.0;
      init_elem.second.upper_bound = 1.0;
      bv_id_interval_stack.push(std::move(init_elem));
    }

    // Processing loop
    OBB<S> bvh_node_bv_with_tf;
    while (!bv_id_interval_stack.empty()) {
      // Current task
      int bv_id = -1;
      Interval<S> bv_parent_interval;
      {
        const auto& stack_top = bv_id_interval_stack.top();
        bv_id = stack_top.first;
        bv_parent_interval = stack_top.second;
        bv_id_interval_stack.pop();
      }

      // Check bv for both leaf and internal
      assert(bv_id < bvh_2->getNumBVs());
      const BVNode<OBB<S>>& bv2_node = bvh_2->getBV(bv_id);

      // Get node2 bv and transform with tf
      const OBB<S>& node2_bv = bv2_node.bv;
      convertBV(node2_bv, tf2, bvh_node_bv_with_tf);

      // OBB ccd test
      Interval<S> obb_interval;
      const bool is_disjoint = BoxPairTranslationalCCD<S>::IsDisjoint(
          shape_bv, shape_bv_displacement, bvh_node_bv_with_tf,
          bv_parent_interval, obb_interval, request.zero_movement_tolerance);
      if (is_disjoint) {
        continue;
      }

      // Not leaf, continue down
      if (!bv2_node.isLeaf()) {
        const int child_1 = bv2_node.leftChild();
        const int child_2 = bv2_node.rightChild();
        bv_id_interval_stack.push({child_1, obb_interval});
        bv_id_interval_stack.push({child_2, obb_interval});
        continue;
      }

      // Leaf case
      assert(bv2_node.isLeaf());
      const int primitive_id = bv2_node.primitiveId();
      const Simplex<S>& simplex = bvh_2->getSimplex(primitive_id);

      // Make the shape solver functor
      ContinuousContactMeta<S> contact;
      contact.o1 = s1;
      contact.o2 = bvh_2;
      contact.b1 = ContinuousContactMeta<S>::kNone;
      contact.b2 = primitive_id;
      contact.external_box_toc = obb_interval;

      // Invoke caller
      ShapePairTranslationalCollisionSolver<S>::template RunShapeSimplex<Shape>(
          *s1, tf1, s1_displacement, simplex, tf2, contact, request, result);

      // Check termination
      if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
        return;
      }
    }
  }

  static void RunExplicitBV(const Shape* s1, const Transform3<S>& tf1,
                            const TranslationalDisplacement<S>& s1_displacement,
                            const BVHModel<OBB<S>>* bvh_2,
                            const Transform3<S>& tf2,
                            const ContinuousCollisionRequest<S>& request,
                            ContinuousCollisionResult<S>& result) {
    if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
      return;
    }

    // Compute bv and displacement in bv
    OBB<S> shape_swept_bv;
    {
      computeBV(*s1, tf1, shape_swept_bv);
      OBB<S> shape_bv_after_displacement = shape_swept_bv;
      shape_bv_after_displacement.To +=
          s1_displacement.scalar_displacement *
          (tf1.linear() * s1_displacement.unit_axis_in_shape1);
      shape_swept_bv += shape_bv_after_displacement;
    }

    // The task stack contains only the bv id
    OBB<S> bvh_node_bv_with_tf;
    std::stack<int> bv_id_stack;
    bv_id_stack.push(0);
    while (!bv_id_stack.empty()) {
      // Current task
      const int bv_id = bv_id_stack.top();
      bv_id_stack.pop();

      // Check bv for both leaf and internal
      assert(bv_id < bvh_2->getNumBVs());
      const BVNode<OBB<S>>& bv2_node = bvh_2->getBV(bv_id);

      // Get node2 bv and transform with tf
      const OBB<S>& node2_bv = bv2_node.bv;
      convertBV(node2_bv, tf2, bvh_node_bv_with_tf);

      // OBB ccd test
      if (!shape_swept_bv.overlap(bvh_node_bv_with_tf)) {
        continue;
      }

      // Not leaf, continue down
      if (!bv2_node.isLeaf()) {
        const int child_1 = bv2_node.leftChild();
        const int child_2 = bv2_node.rightChild();
        bv_id_stack.push(child_1);
        bv_id_stack.push(child_2);
        continue;
      }

      // Leaf case
      assert(bv2_node.isLeaf());
      const int primitive_id = bv2_node.primitiveId();
      const Simplex<S>& simplex = bvh_2->getSimplex(primitive_id);

      // Make the shape solver functor
      ContinuousContactMeta<S> contact;
      contact.o1 = s1;
      contact.o2 = bvh_2;
      contact.b1 = ContinuousContactMeta<S>::kNone;
      contact.b2 = primitive_id;

      // Invoke caller
      ShapePairTranslationalCollisionSolver<S>::template RunShapeSimplex<Shape>(
          *s1, tf1, s1_displacement, simplex, tf2, contact, request, result);

      // Check termination
      if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
        return;
      }
    }
  }
};

/// For mesh pair
template <typename S, typename BV>
struct TranslationalDisplacementBVH_PairSolverImpl {
  // Does not support general bv
  static void Run(const BVHModel<BV>*, const Transform3<S>&,
                  const TranslationalDisplacement<S>&, const BVHModel<BV>*,
                  const Transform3<S>&, const ContinuousCollisionRequest<S>&,
                  ContinuousCollisionResult<S>&) {
    static_assert(std::is_same<typename BV::S, S>::value,
                  "scalar type must match");
    std::cerr << "Not supported" << std::endl;
  }
};

template <typename S>
struct TranslationalDisplacementBVH_PairSolverImpl<S, AABB<S>> {
  // Does not support general bv
  static void Run(const BVHModel<AABB<S>>* bvh1, const Transform3<S>& tf1,
                  const TranslationalDisplacement<S>& bvh1_displacement,
                  const BVHModel<AABB<S>>* bvh2, const Transform3<S>& tf2,
                  const ContinuousCollisionRequest<S>& request,
                  ContinuousCollisionResult<S>& result) {
    // Check initial isOccupied
    if ((!bvh1->isOccupied()) || (!bvh2->isOccupied())) {
      return;
    }

    // Init the disjoint
    FixedOrientationBoxPairTranslationalCCD<S> disjoint;
    disjoint.Initialize(tf1, tf2, bvh1_displacement);

    // The task stack contains the bv id for both side and checked interval
    struct TaskStackElement {
      int bv_id_1{0};
      int bv_id_2{0};
      Interval<S> interval{};

      // Constructor
      TaskStackElement(int bv_id1, int bv_id2, Interval<S> interval_in)
          : bv_id_1(bv_id1),
            bv_id_2(bv_id2),
            interval(std::move(interval_in)) {}
    };
    std::stack<TaskStackElement> task_stack;
    {
      Interval<S> init_interval;
      init_interval.lower_bound = 0.0;
      init_interval.upper_bound = 1.0;
      task_stack.push(TaskStackElement(0, 0, init_interval));
    }

    // Processing loop
    while (!task_stack.empty()) {
      // Gather the id
      int bv_id_1 = -1;
      int bv_id_2 = -1;
      Interval<S> bv_parent_interval;
      {
        const auto& this_task = task_stack.top();
        bv_id_1 = this_task.bv_id_1;
        bv_id_2 = this_task.bv_id_2;
        bv_parent_interval = this_task.interval;
        task_stack.pop();
      }

      // Gather the node
      assert(bv_id_1 < bvh1->getNumBVs() && bv_id_1 >= 0);
      assert(bv_id_2 < bvh2->getNumBVs() && bv_id_2 >= 0);
      const BVNode<AABB<S>>& bv_node1 = bvh1->getBV(bv_id_1);
      const BVNode<AABB<S>>& bv_node2 = bvh2->getBV(bv_id_2);
      const AABB<S>& bv1 = bv_node1.bv;
      const AABB<S>& bv2 = bv_node2.bv;

      // Check bv
      Interval<S> obb_interval;
      const bool is_disjoint =
          disjoint.IsDisjoint(bv1, bv2, bv_parent_interval, obb_interval,
                              request.zero_movement_tolerance);

      // No intersection, move on
      if (is_disjoint) {
        continue;
      }

      // Need to go down
      const bool is_leaf_1 = bv_node1.isLeaf();
      const bool is_leaf_2 = bv_node2.isLeaf();
      if (is_leaf_1 && is_leaf_2) {
        const int primitive_id1 = bv_node1.primitiveId();
        const int primitive_id2 = bv_node2.primitiveId();
        const Simplex<S>& simplex1 = bvh1->getSimplex(primitive_id1);
        const Simplex<S>& simplex2 = bvh2->getSimplex(primitive_id2);

        // Make the shape solver
        ContinuousContactMeta<S> contact;
        contact.o1 = bvh1;
        contact.o2 = bvh2;
        contact.b1 = primitive_id1;
        contact.b2 = primitive_id2;
        contact.external_box_toc = obb_interval;

        // Invoke caller
        ShapePairTranslationalCollisionSolver<S>::RunSimplexPair(
            simplex1, tf1, bvh1_displacement, simplex2, tf2, contact, request,
            result);

        // Check termination
        if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
          return;
        }
      } else {
        // At least one is not leaf, continue
        const bool continue_on_1 =
            (is_leaf_2) || ((!is_leaf_1) && bv1.size() > bv2.size());
        if (continue_on_1) {
          const int child_1 = bv_node1.leftChild();
          const int child_2 = bv_node1.rightChild();
          task_stack.push(TaskStackElement(child_1, bv_id_2, obb_interval));
          task_stack.push(TaskStackElement(child_2, bv_id_2, obb_interval));
        } else {
          const int child_1 = bv_node2.leftChild();
          const int child_2 = bv_node2.rightChild();
          task_stack.push(TaskStackElement(bv_id_1, child_1, obb_interval));
          task_stack.push(TaskStackElement(bv_id_1, child_2, obb_interval));
        }
      }
    }
  }
};

template <typename S>
struct TranslationalDisplacementBVH_PairSolverImpl<S, OBB<S>> {
  // Does not support general bv
  static void Run(const BVHModel<OBB<S>>* bvh1, const Transform3<S>& tf1,
                  const TranslationalDisplacement<S>& bvh1_displacement,
                  const BVHModel<OBB<S>>* bvh2, const Transform3<S>& tf2,
                  const ContinuousCollisionRequest<S>& request,
                  ContinuousCollisionResult<S>& result) {
    // Check initial isOccupied
    if ((!bvh1->isOccupied()) || (!bvh2->isOccupied())) {
      return;
    }

    // The task stack contains the bv id for both side and checked interval
    struct TaskStackElement {
      int bv_id_1{0};
      int bv_id_2{0};
      Interval<S> interval{};

      // Constructor
      TaskStackElement(int bv_id1, int bv_id2, Interval<S> interval_in)
          : bv_id_1(bv_id1),
            bv_id_2(bv_id2),
            interval(std::move(interval_in)) {}
    };
    std::stack<TaskStackElement> task_stack;
    {
      Interval<S> init_interval;
      init_interval.lower_bound = 0.0;
      init_interval.upper_bound = 1.0;
      task_stack.push(TaskStackElement(0, 0, init_interval));
    }

    // For obb checking
    OBB<S> obb1_tf, obb2_tf;
    const Vector3<S> displacement_axis_world =
        tf1.linear() * bvh1_displacement.unit_axis_in_shape1;

    // Processing loop
    while (!task_stack.empty()) {
      // Gather the id
      int bv_id_1 = -1;
      int bv_id_2 = -1;
      Interval<S> bv_parent_interval;
      {
        const auto& this_task = task_stack.top();
        bv_id_1 = this_task.bv_id_1;
        bv_id_2 = this_task.bv_id_2;
        bv_parent_interval = this_task.interval;
        task_stack.pop();
      }

      // Gather the node
      assert(bv_id_1 < bvh1->getNumBVs());
      assert(bv_id_2 < bvh2->getNumBVs());
      const BVNode<OBB<S>>& bv_node1 = bvh1->getBV(bv_id_1);
      const BVNode<OBB<S>>& bv_node2 = bvh2->getBV(bv_id_2);
      const OBB<S>& bv1 = bv_node1.bv;
      const OBB<S>& bv2 = bv_node2.bv;

      // Check bv
      bool is_disjoint = false;
      Interval<S> obb_interval;
      {
        // Compute bv and the displacement
        convertBV(bv1, tf1, obb1_tf);
        convertBV(bv2, tf2, obb2_tf);
        TranslationalDisplacement<S> obb1_displacement;
        obb1_displacement.unit_axis_in_shape1 =
            obb1_tf.axis.transpose() * displacement_axis_world;
        obb1_displacement.scalar_displacement =
            bvh1_displacement.scalar_displacement;

        // Check intersect
        is_disjoint = BoxPairTranslationalCCD<S>::IsDisjoint(
            obb1_tf, obb1_displacement, obb2_tf, bv_parent_interval,
            obb_interval, request.zero_movement_tolerance);
      }

      // No intersection, move on
      if (is_disjoint) {
        continue;
      }

      // Need to go down
      const bool is_leaf_1 = bv_node1.isLeaf();
      const bool is_leaf_2 = bv_node2.isLeaf();
      if (is_leaf_1 && is_leaf_2) {
        const int primitive_id1 = bv_node1.primitiveId();
        const int primitive_id2 = bv_node2.primitiveId();
        const Simplex<S>& simplex1 = bvh1->getSimplex(primitive_id1);
        const Simplex<S>& simplex2 = bvh2->getSimplex(primitive_id2);

        // Make the shape solver
        ContinuousContactMeta<S> contact;
        contact.o1 = bvh1;
        contact.o2 = bvh2;
        contact.b1 = primitive_id1;
        contact.b2 = primitive_id2;
        contact.external_box_toc = obb_interval;

        // Invoke caller
        ShapePairTranslationalCollisionSolver<S>::RunSimplexPair(
            simplex1, tf1, bvh1_displacement, simplex2, tf2, contact, request,
            result);

        // Check termination
        if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
          return;
        }
      } else {
        // At least one is not leaf, continue
        const bool continue_on_1 =
            (is_leaf_2) || ((!is_leaf_1) && bv1.size() > bv2.size());
        if (continue_on_1) {
          const int child_1 = bv_node1.leftChild();
          const int child_2 = bv_node1.rightChild();
          task_stack.push(TaskStackElement(child_1, bv_id_2, obb_interval));
          task_stack.push(TaskStackElement(child_2, bv_id_2, obb_interval));
        } else {
          const int child_1 = bv_node2.leftChild();
          const int child_2 = bv_node2.rightChild();
          task_stack.push(TaskStackElement(bv_id_1, child_1, obb_interval));
          task_stack.push(TaskStackElement(bv_id_1, child_2, obb_interval));
        }
      }
    }
  }
};

/// External interface
template <typename BV>
template <typename Shape>
void TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunShapeMesh(
    const Shape* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement,
    const BVHModel<BV>* bvh2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  TranslationalDisplacementShapeBVHSolverImpl<S, Shape, BV>::Run(
      s1, tf1, s1_displacement, bvh2, tf2, request, result);
}

template <typename BV>
template <typename Shape>
void TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunMeshShape(
    const BVHModel<BV>* bvh1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement, const Shape* s2,
    const Transform3<S>& tf2, const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  TranslationalDisplacement<S> s2_displacement;
  s2_displacement.scalar_displacement = s1_displacement.scalar_displacement;
  s2_displacement.unit_axis_in_shape1 =
      tf2.linear().transpose() *
      (tf1.linear() * (-s1_displacement.unit_axis_in_shape1));
  TranslationalDisplacementShapeBVHSolverImpl<S, Shape, BV>::Run(
      s2, tf2, s2_displacement, bvh1, tf1, request, result);
}

template <typename BV>
void TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunMeshPair(
    const BVHModel<BV>* bvh1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& bvh1_displacement,
    const BVHModel<BV>* bvh2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  TranslationalDisplacementBVH_PairSolverImpl<S, BV>::Run(
      bvh1, tf1, bvh1_displacement, bvh2, tf2, request, result);
}

}  // namespace detail
}  // namespace fcl