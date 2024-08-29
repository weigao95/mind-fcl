#pragma once

namespace fcl {
namespace detail {

template <typename S>
HeightMapCollisionSolver<S>::HeightMapCollisionSolver(
    const GJKSolver<S>* solver)
    : solver(solver) {
  assert(this->solver != nullptr);
}

template <typename S>
void HeightMapCollisionSolver<S>::HeightMapIntersect(
    const HeightMapCollisionGeometry<S>* hm_1,
    const HeightMapCollisionGeometry<S>* hm_2, const Transform3<S>& tf1,
    const Transform3<S>& tf2, const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (hm_1 == nullptr || hm_2 == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Start the processor
  heightMapPairIntersect(*hm_1, *hm_2, tf1, tf2);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::HeightMapShapeIntersect(
    const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Automatic);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::ShapeHeightMapIntersect(
    const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Automatic);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::HeightMapShapeIntersectByBottomLayer(
    const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Flat);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::ShapeHeightMapIntersectByBottomLayer(
    const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Traversal);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::HeightMapShapeIntersectByTraversal(
    const HeightMapCollisionGeometry<S>* heightmap_geometry, const Shape& s,
    const Transform3<S>& tf_hm, const Transform3<S>& tf_shape,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Flat);
}

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::ShapeHeightMapIntersectByTraversal(
    const Shape& s, const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_shape, const Transform3<S>& tf_hm,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr) return;
  request = &request_in;
  result = &result_in;

  // Directly forward to impl
  heightMapShapeIntersectImpl(*heightmap_geometry, s, tf_hm, tf_shape,
                              ShapeHeightMapCollisionType::Traversal);
}

template <typename S>
void HeightMapCollisionSolver<S>::HeightMapOctreeIntersect(
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Octree2CollisionGeometry<S>* octree, const Transform3<S>& tf_hm,
    const Transform3<S>& tf_octree, const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr || octree == nullptr) return;
  request = &request_in;
  result = &result_in;

  heightMapOctreeIntersect(heightmap_geometry, octree, tf_hm, tf_octree);
}

template <typename S>
void HeightMapCollisionSolver<S>::OctreeHeightMapIntersect(
    const Octree2CollisionGeometry<S>* octree,
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_hm,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr || octree == nullptr) return;
  request = &request_in;
  result = &result_in;

  heightMapOctreeIntersect(heightmap_geometry, octree, tf_hm, tf_octree);
}

template <typename S>
template <typename BV>
void HeightMapCollisionSolver<S>::HeightMapBVHIntersect(
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const BVHModel<BV>* bvh, const Transform3<S>& tf_hm,
    const Transform3<S>& tf_bvh, const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr || bvh == nullptr) return;
  request = &request_in;
  result = &result_in;

  heightMapBVHIntersect(heightmap_geometry, bvh, tf_hm, tf_bvh);
}

template <typename S>
template <typename BV>
void HeightMapCollisionSolver<S>::BVH_HeightMapIntersect(
    const BVHModel<BV>* bvh,
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_bvh, const Transform3<S>& tf_hm,
    const CollisionRequest<S>& request_in,
    CollisionResult<S>& result_in) const {
  if (heightmap_geometry == nullptr || bvh == nullptr) return;
  request = &request_in;
  result = &result_in;

  heightMapBVHIntersect(heightmap_geometry, bvh, tf_hm, tf_bvh);
}

}  // namespace detail
}  // namespace fcl