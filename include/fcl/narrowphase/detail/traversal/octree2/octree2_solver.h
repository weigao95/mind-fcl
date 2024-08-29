//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/utility.h"
#include "fcl/math/bv/utility.h"
#include "fcl/math/fixed_rotation_obb_disjoint.h"
#include "fcl/narrowphase/collision_request.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/shape_pair_intersect.h"

namespace fcl {
namespace detail {

template <typename S>
class CollisionSolverOctree2 {
 private:
  using Octree = octree2::Octree<S>;

  // Data
  const GJKSolver<S>* solver;
  mutable const CollisionRequest<S>* request;
  mutable CollisionResult<S>* result;

 public:
  explicit CollisionSolverOctree2(const GJKSolver<S>* solver_in);

  /// Octree vs Shape
  template <typename Shape>
  void OctreeShapeIntersect(const Octree2CollisionGeometry<S>* octree,
                            const Shape& shape, const Transform3<S>& tf_octree,
                            const Transform3<S>& tf_shape,
                            const CollisionRequest<S>& request_in,
                            CollisionResult<S>& result_in) const;
  template <typename Shape>
  void ShapeOctreeIntersect(const Shape& shape,
                            const Octree2CollisionGeometry<S>* octree,
                            const Transform3<S>& tf_shape,
                            const Transform3<S>& tf_octree,
                            const CollisionRequest<S>& request_in,
                            CollisionResult<S>& result_in) const;
  static std::int64_t encodeOctree2Node(
      std::uint32_t octree_node_vector_index, bool is_leaf_node,
      std::uint8_t inner_child_idx = octree2::kInvalidChildIndex);

 private:
  struct OctreeLeafComputeCache {
    // Already init
    ShapePairIntersectSolver<S> shape_solver;

    // Assign before usage
    OBB<S> obb1, obb2;
    AABB<S> aabb1, aabb2;
    Transform3<S> shape1_tf, shape2_tf;
    Box<S> box1, box2;
    ContactMeta<S> contact_meta;
  };

  template <typename Shape>
  void octreeShapeIntersectImpl(const Octree2CollisionGeometry<S>& octree,
                                const Shape& s, const Transform3<S>& tf_octree,
                                const Transform3<S>& tf_shape,
                                OctreeLeafComputeCache& cache) const;
  template <typename Shape>
  void boxToShapeProcessLeafPair(const Octree2CollisionGeometry<S>& octree,
                                 const Transform3<S>& tf_octree,
                                 const Shape& shape,
                                 const Transform3<S>& tf_shape,
                                 std::int64_t encoded_octree_node_idx,
                                 const AABB<S>& voxel_aabb,
                                 OctreeLeafComputeCache& cache) const;

  /// Collision between heightmap and bvh
 public:
  template <typename BV>
  void OctreeBVHIntersect(const Octree2CollisionGeometry<S>* octree_geometry,
                          const BVHModel<BV>* bvh,
                          const Transform3<S>& tf_octree,
                          const Transform3<S>& tf_bvh,
                          const CollisionRequest<S>& request_in,
                          CollisionResult<S>& result_in) const;
  template <typename BV>
  void BVH_OctreeIntersect(const BVHModel<BV>* bvh,
                           const Octree2CollisionGeometry<S>* octree_geometry,
                           const Transform3<S>& tf_bvh,
                           const Transform3<S>& tf_octree,
                           const CollisionRequest<S>& request_in,
                           CollisionResult<S>& result_in) const;

 private:
  template <typename BV>
  void octreeBVHIntersect(const Octree2CollisionGeometry<S>& octree,
                          const BVHModel<BV>* bvh,
                          const Transform3<S>& tf_octree,
                          const Transform3<S>& tf_bvh,
                          OctreeLeafComputeCache& cache) const;
  void boxToSimplexProcessLeafPair(const Octree2CollisionGeometry<S>& octree,
                                   const Transform3<S>& tf_octree,
                                   const AABB<S>& voxel_aabb,
                                   std::int64_t encoded_octree_node_idx,
                                   const CollisionGeometry<S>* bvh,
                                   const Transform3<S>& tf_bvh,
                                   const Simplex<S>& simplex, intptr_t index_2,
                                   OctreeLeafComputeCache& cache) const;

  /// Collision checking between octree pair
 public:
  void OctreePairIntersect(const Octree2CollisionGeometry<S>* tree1,
                           const Octree2CollisionGeometry<S>* tree2,
                           const Transform3<S>& tf_tree1,
                           const Transform3<S>& tf_tree2,
                           const CollisionRequest<S>& request_in,
                           CollisionResult<S>& result_in);

 private:
  void octreePairIntersect(const Octree2CollisionGeometry<S>& tree1,
                           const Octree2CollisionGeometry<S>& tree2,
                           const Transform3<S>& tf_tree1,
                           const Transform3<S>& tf_tree2,
                           const FixedRotationBoxDisjoint<S>& disjoint,
                           OctreeLeafComputeCache& cache) const;
  bool octreePairTwoLeafNode(
      const Octree2CollisionGeometry<S>& tree1,
      const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
      const Transform3<S>& tf_tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const octree2::OctreeLeafNode& tree1_node,
      const octree2::OctreeLeafNode& tree2_node,
      const FixedRotationBoxDisjoint<S>& disjoint,
      OctreeLeafComputeCache& cache) const;
  bool octreePairInnerNodeWithLeafNode(
      const Octree2CollisionGeometry<S>& tree1,
      const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
      const Transform3<S>& tf_tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const octree2::OctreeLeafNode& tree2_node, bool reverse_tree12,
      const FixedRotationBoxDisjoint<S>& disjoint,
      OctreeLeafComputeCache& cache) const;
  void octreePairInnerNodePairAsLeaf(
      const Octree2CollisionGeometry<S>& tree1,
      const Octree2CollisionGeometry<S>& tree2, const Transform3<S>& tf_tree1,
      const Transform3<S>& tf_tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const FixedRotationBoxDisjoint<S>& disjoint,
      OctreeLeafComputeCache& cache) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver-inl.h"
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver_leaf-inl.h"
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver_traverse-inl.h"
