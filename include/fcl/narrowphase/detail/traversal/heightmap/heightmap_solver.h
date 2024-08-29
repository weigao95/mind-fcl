#pragma once

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/heightmap/heightmap_collision_geometry.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/math/fixed_rotation_obb_disjoint.h"
#include "fcl/narrowphase/collision_request.h"
#include "fcl/narrowphase/detail/gjk_solver.h"

namespace fcl {
namespace detail {

/// HeightMapCollisionSolver handle the collision detection between
/// a heightmap with a shape/octree/bvh/heightmap. The algorithm is
/// equivalent to perform collision detection of each bin (box) in
/// the heightmap with another geometry.
///
/// Traversal Algorithm:
/// Except the flatHeightMapShapeIntersectImpl method, all other ones
/// are built upon the traversal of the "HeightMapTree" defined
/// implicitly in the layered heightmap as:
///    1. Each (layer, pixel) pair is a tree node
///     1.1 Each node has 4 children in the next layer (layer + 1)
///     1.2 In a 2x2 expansion scheme, the pixel of the child is
///         (2 * pixel.x + 0/1, 2 * pixel.y + 0/1)
///    2. The root a the global bounding box (not explicitly stored
///       as a layer). Alternatively, we can think of all nodes in
///       the top layer are roots in a "HeightMapForest"
/// Each node in the "HeightMapTree" is a pixel in its layer, thus
/// it represents a box. This box is an AABB for all its children.
///
/// In this way, the traversal algorithm for LayeredHeightMap is
/// very similar to the one in Octree.
///
/// Flat Algorithm:
/// In flatHeightMapShapeIntersectImpl, an heightmap-frame AABB of
/// the geometry is first computed. Then, xOy component of the shape
/// AABB is intersect with the heightmap xOy range to compute a
/// region of interest. All heightmap pixel (box) in that region is
/// collision checked with the given geometry.
/// This method would be suitable when the region of interest is small.
template <typename S>
class HeightMapCollisionSolver {
 private:
  const GJKSolver<S>* solver;
  mutable const CollisionRequest<S>* request;
  mutable CollisionResult<S>* result;
  using FlatHeightMap = heightmap::FlatHeightMap<S>;
  using LayeredHeightMap = heightmap::LayeredHeightMap<S>;

  // Internal option code
  enum class ShapeHeightMapCollisionType { Flat, Traversal, Automatic };

 public:
  explicit HeightMapCollisionSolver(const GJKSolver<S>* solver);

  /// Collision between two heightmaps
  void HeightMapIntersect(const HeightMapCollisionGeometry<S>* hm_1,
                          const HeightMapCollisionGeometry<S>* hm_2,
                          const Transform3<S>& tf1, const Transform3<S>& tf2,
                          const CollisionRequest<S>& request_in,
                          CollisionResult<S>& result_in) const;

  /// Collision between heightmap and shape
  template <typename Shape>
  void HeightMapShapeIntersect(
      const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
      const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  template <typename Shape>
  void ShapeHeightMapIntersect(
      const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  // Specialized algorithm
  template <typename Shape>
  void HeightMapShapeIntersectByBottomLayer(
      const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
      const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  template <typename Shape>
  void ShapeHeightMapIntersectByBottomLayer(
      const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  template <typename Shape>
  void HeightMapShapeIntersectByTraversal(
      const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
      const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  template <typename Shape>
  void ShapeHeightMapIntersectByTraversal(
      const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;

  /// Collision between heightmap and octree
  void HeightMapOctreeIntersect(
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_hm,
      const Transform3<S>& tf_octree, const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  void OctreeHeightMapIntersect(
      const Octree2CollisionGeometry<S>* octree,
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_octree, const Transform3<S>& tf_hm,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;

  /// Collision between heightmap and bvh
  template <typename BV>
  void HeightMapBVHIntersect(
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const BVHModel<BV>* bvh, const Transform3<S>& tf_hm,
      const Transform3<S>& tf_bvh, const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;
  template <typename BV>
  void BVH_HeightMapIntersect(
      const BVHModel<BV>* bvh,
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_bvh, const Transform3<S>& tf_hm,
      const CollisionRequest<S>& request_in,
      CollisionResult<S>& result_in) const;

 private:
  // HeightMap vs Shape
  template <typename Shape>
  void heightMapShapeIntersectImpl(
      const HeightMapCollisionGeometry<S>& heightmap, const Shape& s,
      const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
      ShapeHeightMapCollisionType option_code =
          ShapeHeightMapCollisionType::Automatic) const;

  // Actual impl
  template <typename Shape>
  void flatHeightMapShapeIntersectImpl(
      const HeightMapCollisionGeometry<S>& heightmap, const Shape& shape,
      const AABB<S>& shape_aabb_in_hm,
      const heightmap::PixelSpaceROI& bottom_roi, const Transform3<S>& tf_hm,
      const Transform3<S>& tf_shape) const;
  template <typename Shape>
  void traversalHeightMapShapeIntersectImpl(
      const HeightMapCollisionGeometry<S>& heightmap, const Shape& shape,
      const AABB<S>& shape_aabb_in_hm,
      const heightmap::PixelSpaceROI& bottom_roi, const Transform3<S>& tf_hm,
      const Transform3<S>& tf_shape) const;

  // HeightMap vs itself, octree and bvh
  void heightMapPairIntersect(
      const HeightMapCollisionGeometry<S>& hm_1_geometry,
      const HeightMapCollisionGeometry<S>& hm_2_geometry,
      const Transform3<S>& tf1, const Transform3<S>& tf2) const;
  void heightMapOctreeIntersect(
      const HeightMapCollisionGeometry<S>* hm_geometry,
      const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_hm,
      const Transform3<S>& tf_octree) const;
  template <typename BV>
  void heightMapBVHIntersect(const HeightMapCollisionGeometry<S>* hm_geometry,
                             const BVHModel<BV>* bvh,
                             const Transform3<S>& tf_hm,
                             const Transform3<S>& tf_bvh) const;

  // Processing of leaf case
  void boxToBoxProcessLeafPair(
      const fcl::CollisionGeometry<S>* geom_1, const Transform3<S>& tf1,
      const AABB<S>& node_aabb_1, intptr_t index_1,
      const fcl::CollisionGeometry<S>* geom_2, const Transform3<S>& tf2,
      const AABB<S>& node_aabb_2, intptr_t index_2,
      const FixedRotationBoxDisjoint<S>& obb_disjoint) const;
  template <typename Shape>
  void boxToShapeProcessLeafPair(
      const HeightMapCollisionGeometry<S>& heightmap_geometry,
      const Transform3<S>& tf_heightmap, const Shape& s,
      const Transform3<S>& tf_shape, const heightmap::Pixel& pixel,
      const AABB<S>& pixel_aabb) const;
  void boxToTriangleProcessLeafPair(
      const HeightMapCollisionGeometry<S>* heightmap_geometry,
      const Transform3<S>& tf_hm, const AABB<S>& node_aabb_hm, intptr_t index_1,
      const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
      const Simplex<S>& simplex, intptr_t index_2) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver-inl.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver_leaf-inl.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver_traverse-inl.h"