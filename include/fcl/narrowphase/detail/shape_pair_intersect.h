//
// Created by wei on 23-10-13.
//

#pragma once

#include "fcl/geometry/shape/utility.h"
#include "fcl/math/mesh_simplex.h"
#include "fcl/narrowphase/collision_request.h"
#include "fcl/narrowphase/contact_point.h"
#include "fcl/narrowphase/detail/gjk_solver.h"

namespace fcl {
namespace detail {

// Meta (geometry, index) information that would be written into the colliding
// contact. It also contains option for reversing the normal and object o1/o2
template <typename S>
struct ContactMeta {
  // These are from Contact
  const CollisionGeometry<S>* o1{nullptr};
  const CollisionGeometry<S>* o2{nullptr};
  static constexpr intptr_t kNone = -1;
  intptr_t b1{kNone};
  intptr_t b2{kNone};
  AABB<S> o1_bv{};
  AABB<S> o2_bv{};

  // These are options that would not be written into contact
  bool reverse_normal{false};
  bool reverse_o1_and_o2{false};
  void reset();

  // Write to contact
  void writeToContact(const ContactPoint<S>& contact_point,
                      Contact<S>& contact) const;
  void writeToContact(Contact<S>& contact) const;
};

template <typename S_>
struct ShapePairIntersectSolver {
  using S = S_;
  ShapePairIntersectSolver() : gjk_solver(nullptr) {}
  explicit ShapePairIntersectSolver(const GJKSolver<S>* solver)
      : gjk_solver(solver){};

  template <typename Shape1, typename Shape2>
  void ShapeIntersect(const Shape1& s1, const Transform3<S>& tf1,
                      const Shape2& s2, const Transform3<S>& tf2,
                      const CollisionRequest<S>& request,
                      const ContactMeta<S>& contact_meta,
                      CollisionResult<S>& result) const;
  template <typename Shape1, typename Shape2>
  void ShapeIntersect(const Shape1* s1, const Transform3<S>& tf1,
                      const Shape2* s2, const Transform3<S>& tf2,
                      const CollisionRequest<S>& request,
                      CollisionResult<S>& result) const;

  template <typename Shape>
  void ShapeSimplexIntersect(const Shape& s1, const Transform3<S>& tf1,
                             const Simplex<S>& s2, const Transform3<S>& tf2,
                             const CollisionRequest<S>& request,
                             const ContactMeta<S>& contact_meta,
                             CollisionResult<S>& result) const;

  void SimplexIntersect(const Simplex<S>& s1, const Transform3<S>& tf1,
                        const Simplex<S>& s2, const Transform3<S>& tf2,
                        const Matrix3<S>& rotation_2to1,
                        const Vector3<S>& translation_2in1,
                        const CollisionRequest<S>& request,
                        const ContactMeta<S>& contact_meta,
                        CollisionResult<S>& result) const;

  // Use gjk solver internally
  const GJKSolver<S>* gjk_solver;

 private:
  void trianglePairIntersect(const Simplex<S>& s1, const Transform3<S>& tf1,
                             const Simplex<S>& s2, const Transform3<S>& tf2,
                             const Matrix3<S>& rotation_2to1,
                             const Vector3<S>& translation_2in1,
                             const CollisionRequest<S>& request,
                             const ContactMeta<S>& contact_meta_fn,
                             CollisionResult<S>& result) const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/shape_pair_intersect-inl.h"