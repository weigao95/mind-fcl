//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"
#include "fcl/narrowphase/detail/ccd/ccd_request.h"
#include "fcl/narrowphase/detail/ccd/ccd_result.h"
#include "fcl/narrowphase/detail/ccd/gjk_ccd.h"

namespace fcl {
namespace detail {

template <typename S>
struct ContinuousContactMeta {
  // These are from Contact
  const CollisionGeometry<S>* o1{nullptr};
  const CollisionGeometry<S>* o2{nullptr};
  static constexpr intptr_t kNone = -1;
  std::int64_t b1{kNone};
  std::int64_t b2{kNone};
  AABB<S> o1_bv{};
  AABB<S> o2_bv{};

  // These are options that would not be written into contact
  bool reverse_o1_and_o2{false};
  Interval<S> external_box_toc{};
  void reset();

  // Write to contact
  bool external_toc_valid() const;
  void writeToContact(ContinuousCollisionContact<S>& contact) const;
  void writeToContact(ContinuousCollisionContact<S>& contact,
                      const Interval<S>& toc) const;
};

template <typename S_>
struct ShapePairTranslationalCollisionSolver {
  using S = S_;

  // For shape pair
  template <typename Shape1, typename Shape2>
  static void RunShapePair(const Shape1& s1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const Shape2& s2, const Transform3<S>& tf2,
                           const ContinuousContactMeta<S>& contact_meta,
                           const ContinuousCollisionRequest<S>& request,
                           ContinuousCollisionResult<S>& result);
  template <typename Shape1, typename Shape2>
  static void RunShapePair(const Shape1* s1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const Shape2* s2, const Transform3<S>& tf2,
                           const ContinuousCollisionRequest<S>& request,
                           ContinuousCollisionResult<S>& result);

  // For shape simplex
  template <typename Shape>
  static void RunShapeSimplex(
      const Shape& s1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& s1_displacement, const Simplex<S>& s2,
      const Transform3<S>& tf2, const ContinuousContactMeta<S>& contact_meta,
      const ContinuousCollisionRequest<S>& request,
      ContinuousCollisionResult<S>& result);

  // Simplex pair
  static void RunSimplexPair(
      const Simplex<S>& s1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& s1_displacement, const Simplex<S>& s2,
      const Transform3<S>& tf2, const ContinuousContactMeta<S>& contact_meta,
      const ContinuousCollisionRequest<S>& request,
      ContinuousCollisionResult<S>& result);
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/shape_pair_ccd-inl.h"