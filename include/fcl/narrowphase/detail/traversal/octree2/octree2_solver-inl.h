//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
CollisionSolverOctree2<S>::CollisionSolverOctree2(const GJKSolver<S>* solver_in)
    : solver(solver_in) {
  assert(solver != nullptr);
}

template <typename S>
template <typename Shape>
void CollisionSolverOctree2<S>::OctreeShapeIntersect(
    const Octree2CollisionGeometry<S>* octree, const Shape& s,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_shape,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (octree == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  OctreeLeafComputeCache cache;
  cache.shape_solver = ShapePairIntersectSolver<S>(solver);
  octreeShapeIntersectImpl(*octree, s, tf_octree, tf_shape, cache);
}

template <typename S>
template <typename Shape>
void CollisionSolverOctree2<S>::ShapeOctreeIntersect(
    const Shape& s, const Octree2CollisionGeometry<S>* octree,
    const Transform3<S>& tf_shape, const Transform3<S>& tf_octree,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  static_assert(std::is_same<S, typename Shape::S>::value, "Scalar mismatch");
  if (octree == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  OctreeLeafComputeCache cache;
  cache.shape_solver = ShapePairIntersectSolver<S>(solver);
  octreeShapeIntersectImpl(*octree, s, tf_octree, tf_shape, cache);
}

template <typename S>
template <typename BV>
void CollisionSolverOctree2<S>::OctreeBVHIntersect(
    const Octree2CollisionGeometry<S>* octree_geometry, const BVHModel<BV>* bvh,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_bvh,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  static_assert(std::is_same<S, typename BV::S>::value, "Scalar mismatch");
  if (octree_geometry == nullptr || bvh == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Into impl
  OctreeLeafComputeCache cache;
  cache.shape_solver = ShapePairIntersectSolver<S>(solver);
  octreeBVHIntersect(*octree_geometry, bvh, tf_octree, tf_bvh, cache);
}

template <typename S>
template <typename BV>
void CollisionSolverOctree2<S>::BVH_OctreeIntersect(
    const BVHModel<BV>* bvh, const Octree2CollisionGeometry<S>* octree_geometry,
    const Transform3<S>& tf_bvh, const Transform3<S>& tf_octree,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  static_assert(std::is_same<S, typename BV::S>::value, "Scalar mismatch");
  if (octree_geometry == nullptr || bvh == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Into impl
  OctreeLeafComputeCache cache;
  cache.shape_solver = ShapePairIntersectSolver<S>(solver);
  octreeBVHIntersect(*octree_geometry, bvh, tf_octree, tf_bvh, cache);
}

template <typename S>
void CollisionSolverOctree2<S>::OctreePairIntersect(
    const Octree2CollisionGeometry<S>* tree1,
    const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf_tree1,
    const Transform3<S>& tf_tree2, const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) {
  if (tree1 == nullptr || tree2 == nullptr) return;
  request = &request_in;
  result = &result_in;

  FixedRotationBoxDisjoint<S> disjoint;
  disjoint.initialize(tf_tree1, tf_tree2);
  OctreeLeafComputeCache cache;
  cache.shape_solver = ShapePairIntersectSolver<S>(solver);
  octreePairIntersect(*tree1, *tree2, tf_tree1, tf_tree2, disjoint, cache);
}

}  // namespace detail
}  // namespace fcl
