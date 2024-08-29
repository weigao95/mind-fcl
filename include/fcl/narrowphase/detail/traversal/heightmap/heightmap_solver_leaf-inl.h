#pragma once

#include "fcl/narrowphase/detail/shape_pair_intersect.h"

namespace fcl {
namespace detail {

template <typename S>
template <typename Shape>
void HeightMapCollisionSolver<S>::boxToShapeProcessLeafPair(
    const HeightMapCollisionGeometry<S>& heightmap_geometry,
    const Transform3<S>& tf_heightmap, const Shape& s,
    const Transform3<S>& tf_shape, const heightmap::Pixel& pixel,
    const AABB<S>& pixel_aabb) const {
  // Make the box and do collision
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(pixel_aabb, tf_heightmap, box, box_tf);

  // Make the meta and shape solver
  ContactMeta<S> contact;
  contact.o1 = &heightmap_geometry;
  contact.o2 = &s;
  contact.b1 = encodePixel(pixel);
  contact.b2 = Contact<S>::NONE;
  contact.o1_bv = pixel_aabb;
  ShapePairIntersectSolver<S> shape_solver(solver);
  shape_solver.template ShapeIntersect<Box<S>, Shape>(
      box, box_tf, s, tf_shape, *request, contact, *result);
}

template <typename S>
void HeightMapCollisionSolver<S>::boxToBoxProcessLeafPair(
    const fcl::CollisionGeometry<S>* geom_1, const Transform3<S>& tf1,
    const AABB<S>& aabb_1, fcl::intptr_t index_1,
    const fcl::CollisionGeometry<S>* geom_2, const Transform3<S>& tf2,
    const AABB<S>& aabb_2, fcl::intptr_t index_2,
    const FixedRotationBoxDisjoint<S>& obb_disjoint) const {
  ContactMeta<S> contact;
  contact.o1 = geom_1;
  contact.o2 = geom_2;
  contact.b1 = index_1;
  contact.b2 = index_2;
  contact.o1_bv = aabb_1;
  contact.o2_bv = aabb_2;

  if (!request->isPenetrationEnabled()) {
    // Run detection with given obb
    if (!obb_disjoint.isDisjoint(aabb_1, aabb_2, true)) {
      if (result->numContacts() < request->maxNumContacts()) {
        Contact<S> this_contact;
        contact.writeToContact(this_contact);
        result->addContact(this_contact);
      }
    }
  } else {
    Box<S> box1, box2;
    Transform3<S> box1_tf, box2_tf;
    constructBox(aabb_1, tf1, box1, box1_tf);
    constructBox(aabb_2, tf2, box2, box2_tf);

    // Make the shape solver
    ShapePairIntersectSolver<S> shape_solver(solver);
    shape_solver.template ShapeIntersect<Box<S>, Box<S>>(
        box1, box1_tf, box2, box2_tf, *request, contact, *result);
  }
}

template <typename S>
void HeightMapCollisionSolver<S>::boxToTriangleProcessLeafPair(
    const HeightMapCollisionGeometry<S>* heightmap_geometry,
    const Transform3<S>& tf_hm, const AABB<S>& node_aabb_hm, intptr_t index_1,
    const CollisionGeometry<S>* bvh, const Transform3<S>& tf_bvh,
    const Simplex<S>& simplex, intptr_t index_2) const {
  Box<S> box;
  Transform3<S> box_tf;
  constructBox(node_aabb_hm, tf_hm, box, box_tf);

  // Make the meta and shape solver
  ContactMeta<S> contact;
  contact.o1 = heightmap_geometry;
  contact.o2 = bvh;
  contact.b1 = index_1;
  contact.b2 = index_2;
  contact.o1_bv = node_aabb_hm;
  ShapePairIntersectSolver<S> shape_solver(solver);
  shape_solver.template ShapeSimplexIntersect<Box<S>>(
      box, box_tf, simplex, tf_bvh, *request, contact, *result);
}

}  // namespace detail
}  // namespace fcl