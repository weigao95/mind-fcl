//
// Created by Wei Gao on 2024/6/17.
//

#pragma once

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/heightmap/heightmap_collision_geometry.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"
#include "fcl/narrowphase/detail/ccd/ccd_request.h"
#include "fcl/narrowphase/detail/ccd/ccd_result.h"
#include "fcl/narrowphase/detail/ccd/ccd_solver_utility.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"

namespace fcl {
namespace detail {

template <typename S>
class TranslationalDisplacementHeightMapSolver {
 public:
  using FlatHeightMap = heightmap::FlatHeightMap<S>;
  using LayeredHeightMap = heightmap::LayeredHeightMap<S>;

  template <typename Shape>
  void RunShapeHeightMap(const Shape* s1, const Transform3<S>& tf1,
                         const TranslationalDisplacement<S>& s1_displacement,
                         const HeightMapCollisionGeometry<S>* hm2,
                         const Transform3<S>& tf2,
                         const ContinuousCollisionRequest<S>& request,
                         ContinuousCollisionResult<S>& result);
  template <typename Shape>
  void RunHeightMapShape(const HeightMapCollisionGeometry<S>* hm1,
                         const Transform3<S>& tf1,
                         const TranslationalDisplacement<S>& hm1_displacement,
                         const Shape* s2, const Transform3<S>& tf2,
                         const ContinuousCollisionRequest<S>& request,
                         ContinuousCollisionResult<S>& result);

  void RunHeightMapOctree(const HeightMapCollisionGeometry<S>* hm1,
                          const Transform3<S>& tf1,
                          const TranslationalDisplacement<S>& hm1_displacement,
                          const Octree2CollisionGeometry<S>* octree,
                          const Transform3<S>& tf2,
                          const ContinuousCollisionRequest<S>& request,
                          ContinuousCollisionResult<S>& result);
  void RunOctreeHeightMap(
      const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_octree,
      const TranslationalDisplacement<S>& octree_displacement,
      const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
      const ContinuousCollisionRequest<S>& request,
      ContinuousCollisionResult<S>& result);

  void RunHeightMapObbBVH(const HeightMapCollisionGeometry<S>* hm1,
                          const Transform3<S>& tf1,
                          const TranslationalDisplacement<S>& hm1_displacement,
                          const BVHModel<OBB<S>>* bvh2,
                          const Transform3<S>& tf2,
                          const ContinuousCollisionRequest<S>& request,
                          ContinuousCollisionResult<S>& result);
  void RunObbBVH_HeightMap(
      const BVHModel<OBB<S>>* bvh1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& bvh1_displacement,
      const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
      const ContinuousCollisionRequest<S>& request,
      ContinuousCollisionResult<S>& result);

  void RunHeightMapPair(const HeightMapCollisionGeometry<S>* hm1,
                        const Transform3<S>& tf1,
                        const TranslationalDisplacement<S>& hm1_displacement,
                        const HeightMapCollisionGeometry<S>* hm2,
                        const Transform3<S>& tf2,
                        const ContinuousCollisionRequest<S>& request,
                        ContinuousCollisionResult<S>& result);

 private:
  // Data
  mutable const ContinuousCollisionRequest<S>* request{nullptr};
  mutable ContinuousCollisionResult<S>* result{nullptr};

  // Shape vs octree
  template <typename Shape>
  void runShapeHeightMap(const Shape* s1, const Transform3<S>& tf1,
                         const TranslationalDisplacement<S>& s1_displacement,
                         const HeightMapCollisionGeometry<S>* hm2,
                         const Transform3<S>& tf2) const;
  template <typename Shape>
  void shapeToBoxProcessLeafPair(
      const Shape* s1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& s1_displacement,
      const HeightMapCollisionGeometry<S>* hm2, const Transform3<S>& tf2,
      const heightmap::Pixel& pixel, const AABB<S>& pixel_aabb) const;

  void runHeightMapObbBVH(const HeightMapCollisionGeometry<S>* hm1,
                          const Transform3<S>& tf_hm1,
                          const TranslationalDisplacement<S>& hm1_displacement,
                          const BVHModel<OBB<S>>* bvh2,
                          const Transform3<S>& tf_bvh2) const;
  void boxToSimplexProcessLeafPair(
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_hm,
      const TranslationalDisplacement<S>& hm1_displacement,
      const AABB<S>& node_aabb_hm, std::int64_t index_1,
      const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
      const Simplex<S>& simplex, std::int64_t index_2) const;

  void runHeightMapOctree(const HeightMapCollisionGeometry<S>* hm_geometry,
                          const Transform3<S>& tf_hm,
                          const TranslationalDisplacement<S>& hm_displacement,
                          const Octree2CollisionGeometry<S>* octree,
                          const Transform3<S>& tf_octree) const;
  void runHeightMapPair(const HeightMapCollisionGeometry<S>* hm1,
                        const Transform3<S>& tf1,
                        const TranslationalDisplacement<S>& hm1_displacement,
                        const HeightMapCollisionGeometry<S>* hm2,
                        const Transform3<S>& tf2) const;
  void boxToBoxProcessLeafPair(
      const fcl::CollisionGeometry<S>* geom_1, const AABB<S>& node_aabb_1,
      std::int64_t index_1, const fcl::CollisionGeometry<S>* geom_2,
      const AABB<S>& node_aabb_2, std::int64_t index_2,
      const FixedOrientationBoxPairTranslationalCCD<S>& disjoint) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/heightmap_ccd_solver-inl.h"
