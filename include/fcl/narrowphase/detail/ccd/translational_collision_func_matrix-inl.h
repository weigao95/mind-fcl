#pragma once

#include "fcl/narrowphase/detail/ccd/bvh_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/heightmap_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"

namespace fcl {
namespace detail {

// Shape pair
//==============================================================================
template <typename Shape1, typename Shape2>
void ShapePairTranslationalCollideImpl(
    const CollisionGeometry<typename Shape1::S>* o1,
    const Transform3<typename Shape1::S>& tf1,
    const TranslationalDisplacement<typename Shape1::S>& o1_displacement,
    const CollisionGeometry<typename Shape1::S>* o2,
    const Transform3<typename Shape1::S>& tf2,
    const ContinuousCollisionRequest<typename Shape1::S>& request,
    ContinuousCollisionResult<typename Shape1::S>& result) {
  // Obtain the objects
  using S = typename Shape1::S;
  static_assert(std::is_same<typename Shape1::S, typename Shape2::S>::value,
                "scalar type must match");
  const auto* obj1 = static_cast<const Shape1*>(o1);
  const auto* obj2 = static_cast<const Shape2*>(o2);

  // Forward to interface
  using Solver = ShapePairTranslationalCollisionSolver<S>;
  Solver::template RunShapePair<Shape1, Shape2>(obj1, tf1, o1_displacement,
                                                obj2, tf2, request, result);
}

// Shape with BVH
//==============================================================================
template <typename Shape, typename BV>
void ShapeBVH_TranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  static_assert(std::is_same<S, typename BV::S>::value,
                "scalar type must match");
  const auto* obj1 = static_cast<const Shape*>(o1);
  const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOrientedNodeBVHSolver<BV>;
  Solver::template RunShapeMesh<Shape>(obj1, tf1, o1_displacement, obj2, tf2,
                                       request, result);
}

// BVH with Shape
//==============================================================================
template <typename Shape, typename BV>
void BVH_ShapeTranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  static_assert(std::is_same<S, typename BV::S>::value,
                "scalar type must match");
  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOrientedNodeBVHSolver<BV>;
  Solver::template RunMeshShape<Shape>(obj1, tf1, o1_displacement, obj2, tf2,
                                       request, result);
}

// BVH Pair
//==============================================================================
template <typename BV>
void BVH_PairTranslationalCollideImpl(
    const CollisionGeometry<typename BV::S>* o1,
    const Transform3<typename BV::S>& tf1,
    const TranslationalDisplacement<typename BV::S>& o1_displacement,
    const CollisionGeometry<typename BV::S>* o2,
    const Transform3<typename BV::S>& tf2,
    const ContinuousCollisionRequest<typename BV::S>& request,
    ContinuousCollisionResult<typename BV::S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const BVHModel<BV>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<BV>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOrientedNodeBVHSolver<BV>;
  Solver::RunMeshPair(obj1, tf1, o1_displacement, obj2, tf2, request, result);
}

// Shape with Octree
//==============================================================================
template <typename Shape>
void ShapeOctreeTranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  const auto* obj1 = static_cast<const Shape*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOctreeSolver<S>;
  Solver solver;
  solver.template RunShapeOctree<Shape>(obj1, tf1, o1_displacement, obj2, tf2,
                                        request, result);
}

// Shape with Octree
//==============================================================================
template <typename Shape>
void OctreeShapeTranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOctreeSolver<S>;
  Solver solver;
  solver.template RunOctreeShape<Shape>(obj1, tf1, o1_displacement, obj2, tf2,
                                        request, result);
}

// Octree With BVH
//==============================================================================
template <typename S>
void OctreeObbBVH_TranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<OBB<S>>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOctreeSolver<S>;
  Solver solver;
  solver.RunOctreeObbBVH(obj1, tf1, o1_displacement, obj2, tf2, request,
                         result);
}

template <typename S>
void ObbBVH_OctreeTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const BVHModel<OBB<S>>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOctreeSolver<S>;
  Solver solver;
  solver.RunObbBVH_Octree(obj1, tf1, o1_displacement, obj2, tf2, request,
                          result);
}

// Octree Pair
//==============================================================================
template <typename S>
void OctreePairTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementOctreeSolver<S>;
  Solver solver;
  solver.RunOctreePair(obj1, tf1, o1_displacement, obj2, tf2, request, result);
}

// Shape with HeightMap
//==============================================================================
template <typename Shape>
void ShapeHeightMapTranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  const auto* obj1 = static_cast<const Shape*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.template RunShapeHeightMap<Shape>(obj1, tf1, o1_displacement, obj2,
                                           tf2, request, result);
}

template <typename Shape>
void HeightMapShapeTranslationalCollideImpl(
    const CollisionGeometry<typename Shape::S>* o1,
    const Transform3<typename Shape::S>& tf1,
    const TranslationalDisplacement<typename Shape::S>& o1_displacement,
    const CollisionGeometry<typename Shape::S>* o2,
    const Transform3<typename Shape::S>& tf2,
    const ContinuousCollisionRequest<typename Shape::S>& request,
    ContinuousCollisionResult<typename Shape::S>& result) {
  // Obtain the objects
  using S = typename Shape::S;
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Shape*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.template RunHeightMapShape<Shape>(obj1, tf1, o1_displacement, obj2,
                                           tf2, request, result);
}

// HeightMap BVH
//==============================================================================
template <typename S>
void HeightMapObbBVH_TranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const BVHModel<OBB<S>>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.RunHeightMapObbBVH(obj1, tf1, o1_displacement, obj2, tf2, request,
                            result);
}

template <typename S>
void ObbBVH_HeightMapTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const BVHModel<OBB<S>>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.RunObbBVH_HeightMap(obj1, tf1, o1_displacement, obj2, tf2, request,
                             result);
}

// HeightMap Octree
//==============================================================================
template <typename S>
void HeightMapOctreeTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const Octree2CollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.RunHeightMapOctree(obj1, tf1, o1_displacement, obj2, tf2, request,
                            result);
}

template <typename S>
void OctreeHeightMapTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const Octree2CollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.RunOctreeHeightMap(obj1, tf1, o1_displacement, obj2, tf2, request,
                            result);
}

// HeightMap Pair
//==============================================================================
template <typename S>
void HeightMapPairTranslationalCollideImpl(
    const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& o1_displacement,
    const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Obtain the objects
  const auto* obj1 = static_cast<const HeightMapCollisionGeometry<S>*>(o1);
  const auto* obj2 = static_cast<const HeightMapCollisionGeometry<S>*>(o2);

  // Forward to interface
  using Solver = TranslationalDisplacementHeightMapSolver<S>;
  Solver solver;
  solver.RunHeightMapPair(obj1, tf1, o1_displacement, obj2, tf2, request,
                          result);
}

template <typename S>
TranslationalCollisionFunctionMatrix<
    S>::TranslationalCollisionFunctionMatrix() {
  // Empty init
  for (int i = 0; i < NODE_COUNT; ++i) {
    for (int j = 0; j < NODE_COUNT; ++j) {
      collision_matrix[i][j] = nullptr;
    }
  }

  /// Shape vs shape
  // clang-format off
  collision_matrix[GEOM_BOX][GEOM_BOX] = &ShapePairTranslationalCollideImpl<Box<S>, Box<S>>;
  collision_matrix[GEOM_BOX][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Box<S>, Sphere<S>>;
  collision_matrix[GEOM_BOX][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Box<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_BOX][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Box<S>, Capsule<S>>;
  collision_matrix[GEOM_BOX][GEOM_CONE] = &ShapePairTranslationalCollideImpl<Box<S>, Cone<S>>;
  collision_matrix[GEOM_BOX][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Box<S>, Cylinder<S>>;
  collision_matrix[GEOM_BOX][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Box<S>, Convex<S>>;
  collision_matrix[GEOM_BOX][GEOM_PLANE] = &ShapePairTranslationalCollideImpl<Box<S>, Plane<S>>;
  collision_matrix[GEOM_BOX][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Box<S>, Halfspace<S>>;

  collision_matrix[GEOM_SPHERE][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Sphere<S>, Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Sphere<S>, Sphere<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Sphere<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Sphere<S>, Capsule<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Sphere<S>, Cone<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Sphere<S>, Cylinder<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Sphere<S>, Convex<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Sphere<S>, Plane<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Sphere<S>, Halfspace<S>>;

  collision_matrix[GEOM_ELLIPSOID][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Box<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Capsule<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Cone<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Cylinder<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Convex<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Plane<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Ellipsoid<S>, Halfspace<S>>;

  collision_matrix[GEOM_CAPSULE][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Capsule<S>, Box<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Capsule<S>, Sphere<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Capsule<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Capsule<S>, Capsule<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Capsule<S>, Cone<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Capsule<S>, Cylinder<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Capsule<S>, Convex<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Capsule<S>, Plane<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Capsule<S>, Halfspace<S>>;

  collision_matrix[GEOM_CONE][GEOM_BOX] = &ShapePairTranslationalCollideImpl<Cone<S>, Box<S>>;
  collision_matrix[GEOM_CONE][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Cone<S>, Sphere<S>>;
  collision_matrix[GEOM_CONE][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Cone<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CONE][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Cone<S>, Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_CONE] = &ShapePairTranslationalCollideImpl<Cone<S>, Cone<S>>;
  collision_matrix[GEOM_CONE][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Cone<S>, Cylinder<S>>;
  collision_matrix[GEOM_CONE][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Cone<S>, Convex<S>>;
  collision_matrix[GEOM_CONE][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Cone<S>, Plane<S>>;
  collision_matrix[GEOM_CONE][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Cone<S>, Halfspace<S>>;

  collision_matrix[GEOM_CYLINDER][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Box<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Sphere<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Capsule<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Cylinder<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Convex<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Plane<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Cylinder<S>, Halfspace<S>>;

  collision_matrix[GEOM_CONVEX][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Convex<S>, Box<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Convex<S>, Sphere<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Convex<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Convex<S>, Capsule<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Convex<S>, Cone<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Convex<S>, Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Convex<S>, Convex<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Convex<S>, Plane<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Convex<S>, Halfspace<S>>;

  collision_matrix[GEOM_PLANE][GEOM_BOX] = &ShapePairTranslationalCollideImpl<Plane<S>, Box<S>>;
  collision_matrix[GEOM_PLANE][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Plane<S>, Sphere<S>>;
  collision_matrix[GEOM_PLANE][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Plane<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Plane<S>, Capsule<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Plane<S>, Cone<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Plane<S>, Cylinder<S>>;
  collision_matrix[GEOM_PLANE][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Plane<S>, Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Plane<S>, Plane<S>>;
  collision_matrix[GEOM_PLANE][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Plane<S>, Halfspace<S>>;

  collision_matrix[GEOM_HALFSPACE][GEOM_BOX] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Box<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_SPHERE] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Sphere<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_ELLIPSOID] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Ellipsoid<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CAPSULE] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Capsule<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CONE] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Cone<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CYLINDER] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Cylinder<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_CONVEX] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Convex<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_PLANE] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_HALFSPACE] =&ShapePairTranslationalCollideImpl<Halfspace<S>, Halfspace<S>>;

  /// BVH
  collision_matrix[GEOM_BOX][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Box<S>, AABB<S>>;
  collision_matrix[GEOM_SPHERE][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Sphere<S>, AABB<S>>;
  collision_matrix[GEOM_ELLIPSOID][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Ellipsoid<S>, AABB<S>>;
  collision_matrix[GEOM_CAPSULE][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Capsule<S>, AABB<S>>;
  collision_matrix[GEOM_CONE][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Cone<S>, AABB<S>>;
  collision_matrix[GEOM_CYLINDER][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Cylinder<S>, AABB<S>>;
  collision_matrix[GEOM_CONVEX][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Convex<S>, AABB<S>>;
  collision_matrix[GEOM_PLANE][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Plane<S>, AABB<S>>;
  collision_matrix[GEOM_HALFSPACE][BV_AABB] = &ShapeBVH_TranslationalCollideImpl<Halfspace<S>, AABB<S>>;

  collision_matrix[BV_AABB][GEOM_BOX] = &BVH_ShapeTranslationalCollideImpl<Box<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_SPHERE] = &BVH_ShapeTranslationalCollideImpl<Sphere<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_ELLIPSOID] = &BVH_ShapeTranslationalCollideImpl<Ellipsoid<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_CAPSULE] = &BVH_ShapeTranslationalCollideImpl<Capsule<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_CONE] = &BVH_ShapeTranslationalCollideImpl<Cone<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_CYLINDER] = &BVH_ShapeTranslationalCollideImpl<Cylinder<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_CONVEX] = &BVH_ShapeTranslationalCollideImpl<Convex<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_PLANE] = &BVH_ShapeTranslationalCollideImpl<Plane<S>, AABB<S>>;
  collision_matrix[BV_AABB][GEOM_HALFSPACE] = &BVH_ShapeTranslationalCollideImpl<Halfspace<S>, AABB<S>>;

  collision_matrix[GEOM_BOX][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Box<S>, OBB<S>>;
  collision_matrix[GEOM_SPHERE][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Sphere<S>, OBB<S>>;
  collision_matrix[GEOM_ELLIPSOID][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Ellipsoid<S>, OBB<S>>;
  collision_matrix[GEOM_CAPSULE][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Capsule<S>, OBB<S>>;
  collision_matrix[GEOM_CONE][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Cone<S>, OBB<S>>;
  collision_matrix[GEOM_CYLINDER][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Cylinder<S>, OBB<S>>;
  collision_matrix[GEOM_CONVEX][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Convex<S>, OBB<S>>;
  collision_matrix[GEOM_PLANE][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Plane<S>, OBB<S>>;
  collision_matrix[GEOM_HALFSPACE][BV_OBB] = &ShapeBVH_TranslationalCollideImpl<Halfspace<S>, OBB<S>>;

  collision_matrix[BV_OBB][GEOM_BOX] = &BVH_ShapeTranslationalCollideImpl<Box<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_SPHERE] = &BVH_ShapeTranslationalCollideImpl<Sphere<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_ELLIPSOID] = &BVH_ShapeTranslationalCollideImpl<Ellipsoid<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_CAPSULE] = &BVH_ShapeTranslationalCollideImpl<Capsule<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_CONE] = &BVH_ShapeTranslationalCollideImpl<Cone<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_CYLINDER] = &BVH_ShapeTranslationalCollideImpl<Cylinder<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_CONVEX] = &BVH_ShapeTranslationalCollideImpl<Convex<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_PLANE] = &BVH_ShapeTranslationalCollideImpl<Plane<S>, OBB<S>>;
  collision_matrix[BV_OBB][GEOM_HALFSPACE] = &BVH_ShapeTranslationalCollideImpl<Halfspace<S>, OBB<S>>;

  collision_matrix[BV_OBB][BV_OBB] = &BVH_PairTranslationalCollideImpl<OBB<S>>;
  collision_matrix[BV_AABB][BV_AABB] = &BVH_PairTranslationalCollideImpl<AABB<S>>;

  /// Octree
  collision_matrix[GEOM_BOX][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_OCTREE2] = &ShapeOctreeTranslationalCollideImpl<Halfspace<S>>;

  collision_matrix[GEOM_OCTREE2][GEOM_BOX] = &OctreeShapeTranslationalCollideImpl<Box<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_SPHERE] =&OctreeShapeTranslationalCollideImpl<Sphere<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_ELLIPSOID] =&OctreeShapeTranslationalCollideImpl<Ellipsoid<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CAPSULE] =&OctreeShapeTranslationalCollideImpl<Capsule<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CONE] =&OctreeShapeTranslationalCollideImpl<Cone<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CYLINDER] =&OctreeShapeTranslationalCollideImpl<Cylinder<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_CONVEX] =&OctreeShapeTranslationalCollideImpl<Convex<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_PLANE] =&OctreeShapeTranslationalCollideImpl<Plane<S>>;
  collision_matrix[GEOM_OCTREE2][GEOM_HALFSPACE] =&OctreeShapeTranslationalCollideImpl<Halfspace<S>>;

  collision_matrix[GEOM_OCTREE2][GEOM_OCTREE2] =&OctreePairTranslationalCollideImpl<S>;

  collision_matrix[GEOM_OCTREE2][BV_OBB] = &OctreeObbBVH_TranslationalCollideImpl<S>;
  collision_matrix[BV_OBB][GEOM_OCTREE2] = &ObbBVH_OctreeTranslationalCollideImpl<S>;

  /// HeightMap
  collision_matrix[GEOM_BOX][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Box<S>>;
  collision_matrix[GEOM_SPHERE][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Sphere<S>>;
  collision_matrix[GEOM_ELLIPSOID][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Ellipsoid<S>>;
  collision_matrix[GEOM_CAPSULE][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Capsule<S>>;
  collision_matrix[GEOM_CONE][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Cone<S>>;
  collision_matrix[GEOM_CYLINDER][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Cylinder<S>>;
  collision_matrix[GEOM_CONVEX][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Convex<S>>;
  collision_matrix[GEOM_PLANE][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Plane<S>>;
  collision_matrix[GEOM_HALFSPACE][GEOM_HEIGHTMAP] = &ShapeHeightMapTranslationalCollideImpl<Halfspace<S>>;

  collision_matrix[GEOM_HEIGHTMAP][GEOM_BOX] = &HeightMapShapeTranslationalCollideImpl<Box<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_SPHERE] =&HeightMapShapeTranslationalCollideImpl<Sphere<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_ELLIPSOID] =&HeightMapShapeTranslationalCollideImpl<Ellipsoid<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CAPSULE] =&HeightMapShapeTranslationalCollideImpl<Capsule<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CONE] =&HeightMapShapeTranslationalCollideImpl<Cone<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CYLINDER] =&HeightMapShapeTranslationalCollideImpl<Cylinder<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_CONVEX] =&HeightMapShapeTranslationalCollideImpl<Convex<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_PLANE] =&HeightMapShapeTranslationalCollideImpl<Plane<S>>;
  collision_matrix[GEOM_HEIGHTMAP][GEOM_HALFSPACE] =&HeightMapShapeTranslationalCollideImpl<Halfspace<S>>;

  collision_matrix[GEOM_HEIGHTMAP][GEOM_OCTREE2] = &HeightMapOctreeTranslationalCollideImpl<S>;
  collision_matrix[GEOM_OCTREE2][GEOM_HEIGHTMAP] = &OctreeHeightMapTranslationalCollideImpl<S>;
  collision_matrix[GEOM_HEIGHTMAP][BV_OBB] = &HeightMapObbBVH_TranslationalCollideImpl<S>;
  collision_matrix[BV_OBB][GEOM_HEIGHTMAP] = &ObbBVH_HeightMapTranslationalCollideImpl<S>;

  collision_matrix[GEOM_HEIGHTMAP][GEOM_HEIGHTMAP] = &HeightMapPairTranslationalCollideImpl<S>;
  // clang-format on
}

}  // namespace detail
}  // namespace fcl
