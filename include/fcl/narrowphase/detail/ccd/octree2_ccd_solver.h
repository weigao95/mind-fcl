//
// Created by wei on 24-6-15.
//

#pragma once

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/geometry/shape/utility.h"
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"
#include "fcl/narrowphase/detail/ccd/ccd_solver_utility.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver.h"

namespace fcl {
namespace detail {

template <typename S>
class TranslationalDisplacementOctreeSolver {
 public:
  // Typedef
  using Octree = octree2::Octree<S>;
  explicit TranslationalDisplacementOctreeSolver() = default;
  ~TranslationalDisplacementOctreeSolver() = default;

  // Shape vs octree
  template <typename Shape>
  void RunShapeOctree(const Shape* s1, const Transform3<S>& tf1,
                      const TranslationalDisplacement<S>& s1_displacement,
                      const Octree2CollisionGeometry<S>* tree2,
                      const Transform3<S>& tf2,
                      const ContinuousCollisionRequest<S>& request,
                      ContinuousCollisionResult<S>& result) const;
  template <typename Shape>
  void RunOctreeShape(const Octree2CollisionGeometry<S>* tree1,
                      const Transform3<S>& tf1,
                      const TranslationalDisplacement<S>& tree1_displacement,
                      const Shape* s2, const Transform3<S>& tf2,
                      const ContinuousCollisionRequest<S>& request,
                      ContinuousCollisionResult<S>& result) const;

  // Octree with bvh
  void RunOctreeObbBVH(const Octree2CollisionGeometry<S>* tree1,
                       const Transform3<S>& tf1,
                       const TranslationalDisplacement<S>& tree1_displacement,
                       const BVHModel<OBB<S>>* tree2, const Transform3<S>& tf2,
                       const ContinuousCollisionRequest<S>& request,
                       ContinuousCollisionResult<S>& result);
  void RunObbBVH_Octree(const BVHModel<OBB<S>>* tree1, const Transform3<S>& tf1,
                        const TranslationalDisplacement<S>& tree1_displacement,
                        const Octree2CollisionGeometry<S>* tree2,
                        const Transform3<S>& tf2,
                        const ContinuousCollisionRequest<S>& request,
                        ContinuousCollisionResult<S>& result);

  // Octree pair
  void RunOctreePair(const Octree2CollisionGeometry<S>* tree1,
                     const Transform3<S>& tf1,
                     const TranslationalDisplacement<S>& tree1_displacement,
                     const Octree2CollisionGeometry<S>* tree2,
                     const Transform3<S>& tf2,
                     const ContinuousCollisionRequest<S>& request,
                     ContinuousCollisionResult<S>& result);

 private:
  // Data
  mutable const ContinuousCollisionRequest<S>* request{nullptr};
  mutable ContinuousCollisionResult<S>* result{nullptr};

  // Shape vs octree
  template <typename Shape>
  void runShapeOctree(const Shape* s1, const Transform3<S>& tf1,
                      const TranslationalDisplacement<S>& s1_displacement,
                      const Octree2CollisionGeometry<S>* tree2,
                      const Transform3<S>& tf2) const;
  template <typename Shape>
  void shapeToBoxProcessLeafPair(
      const Shape* s1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& s1_displacement,
      const Octree2CollisionGeometry<S>* tree2, const Transform3<S>& tf2,
      std::int64_t encoded_octree_node_idx, const AABB<S>& voxel_aabb) const;

  // Octree with bvh
  void runOctreeObbBVH(const Octree2CollisionGeometry<S>* tree1,
                       const Transform3<S>& tf1,
                       const TranslationalDisplacement<S>& s1_displacement,
                       const BVHModel<OBB<S>>* bvh2,
                       const Transform3<S>& tf2) const;
  void boxToSimplexProcessLeafPair(
      const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_octree,
      const TranslationalDisplacement<S>& octree_displacement,
      const AABB<S>& voxel_aabb, std::int64_t encoded_octree_node_idx,
      const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
      const Simplex<S>& simplex, intptr_t index_2) const;

  // Octree pair
  void runOctreePair(
      const Octree2CollisionGeometry<S>* tree1,
      const Octree2CollisionGeometry<S>* tree2,
      const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const;
  bool runOctreePairTwoLeafNode(
      const Octree2CollisionGeometry<S>* tree1,
      const Octree2CollisionGeometry<S>* tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const octree2::OctreeLeafNode& tree1_node,
      const octree2::OctreeLeafNode& tree2_node,
      const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const;
  bool runOctreePairInnerNodeWithLeafNode(
      const Octree2CollisionGeometry<S>* tree1,
      const Octree2CollisionGeometry<S>* tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const octree2::OctreeLeafNode& tree2_node, bool reverse_tree12,
      const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const;
  void runOctreePairInnerNodePairAsLeaf(
      const Octree2CollisionGeometry<S>* tree1,
      const Octree2CollisionGeometry<S>* tree2,
      const octree2::OctreeTraverseStackElement<S>& tree1_elem,
      const octree2::OctreeTraverseStackElement<S>& tree2_elem,
      const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver-inl.h"