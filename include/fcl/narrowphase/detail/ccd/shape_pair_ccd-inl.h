#pragma once

namespace fcl {
namespace detail {

template <typename S>
void ContinuousContactMeta<S>::reset() {
  o1 = nullptr;
  o2 = nullptr;
  b1 = kNone;
  b2 = kNone;
  reverse_o1_and_o2 = false;
  external_box_toc.lower_bound = external_box_toc.upper_bound = S(-1.0);
}

template <typename S>
bool ContinuousContactMeta<S>::external_toc_valid() const {
  return external_box_toc.lower_bound >= 0 && external_box_toc.upper_bound >= 0;
}

template <typename S>
void ContinuousContactMeta<S>::writeToContact(
    ContinuousCollisionContact<S>& contact) const {
  if (!reverse_o1_and_o2) {
    contact.o1 = o1;
    contact.o2 = o2;
    contact.b1 = b1;
    contact.b2 = b2;
    contact.o1_bv = o1_bv;
    contact.o2_bv = o2_bv;
  } else {
    contact.o1 = o2;
    contact.o2 = o1;
    contact.b1 = b2;
    contact.b2 = b1;
    contact.o1_bv = o2_bv;
    contact.o2_bv = o1_bv;
  }

  // Write my toc
  contact.toc = external_box_toc;
}

template <typename S>
void ContinuousContactMeta<S>::writeToContact(
    ContinuousCollisionContact<S>& contact, const Interval<S>& toc) const {
  writeToContact(contact);
  if (toc.lower_bound >= 0 && toc.upper_bound >= 0) {
    contact.toc = toc;
  }
}

template <typename S, typename Shape1, typename Shape2>
struct ShapePairTranslationalCollisionImpl {
  static bool RunIntersect(const Shape1& s1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const Shape2& s2, const Transform3<S>& tf2,
                           const ContinuousCollisionRequest<S>& request,
                           TimeOfCollisionRequestType standalone_request_type,
                           Interval<S>& toc) {
    static_assert(std::is_same<S, typename Shape1::S>::value,
                  "scalar type must match");
    static_assert(std::is_same<S, typename Shape2::S>::value,
                  "scalar type must match");
    // Not occupied
    if ((!s1.isOccupied()) || (!s2.isOccupied())) {
      return false;
    }

    // Make gjk shape
    const auto gjk_s1 = constructGJKGeometry(&s1);
    const auto gjk_s2 = constructGJKGeometry(&s2);

    // Depends on request type
    if (standalone_request_type == TimeOfCollisionRequestType::kNotRequested) {
      // Binary collision
      const auto is_intersect =
          TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionBinary(
              gjk_s1, tf1, s1_displacement, gjk_s2, tf2,
              request.max_gjk_iterations, request.gjk_tolerance);

      // Write the contact
      toc.lower_bound = toc.upper_bound = -1.0;
      return is_intersect;
    } else if (standalone_request_type ==
               TimeOfCollisionRequestType::kOneTocSample) {
      const auto output =
          TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionOneTocSample(
              gjk_s1, tf1, s1_displacement, gjk_s2, tf2,
              request.max_gjk_iterations, request.gjk_tolerance);

      // Write the contact
      toc.lower_bound = toc.upper_bound = output.second;
      return output.first;
    } else {
      assert(standalone_request_type ==
             TimeOfCollisionRequestType::kBoxApproximate);

      // Pre-check with box
      OBB<S> o1_bv;
      OBB<S> o2_bv;
      fcl::convertBV(s1.aabb_local, tf1, o1_bv);
      fcl::convertBV(s2.aabb_local, tf2, o2_bv);
      const auto obb_disjoint = BoxPairTranslationalCCD<S>::IsDisjoint(
          o1_bv, s1_displacement, o2_bv, toc, request.zero_movement_tolerance);
      if (obb_disjoint) return false;

      // Box checking passed, now use gjk
      return TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionBinary(
          gjk_s1, tf1, s1_displacement, gjk_s2, tf2, request.max_gjk_iterations,
          request.gjk_tolerance);
    }
  }
};

template <typename S>
struct ShapePairTranslationalCollisionImpl<S, Box<S>, Box<S>> {
  static bool RunIntersect(const Box<S>& s1, const Transform3<S>& tf1,
                           const TranslationalDisplacement<S>& s1_displacement,
                           const Box<S>& s2, const Transform3<S>& tf2,
                           const ContinuousCollisionRequest<S>& request,
                           TimeOfCollisionRequestType, Interval<S>& toc) {
    // Debug Check
    assert(s1.aabb_local.max_[0] >= s1.aabb_local.min_[0]);
    assert(s1.aabb_local.max_[1] >= s1.aabb_local.min_[1]);
    assert(s1.aabb_local.max_[2] >= s1.aabb_local.min_[2]);

    // Make the obb
    OBB<S> o1_bv;
    OBB<S> o2_bv;
    computeBV(s1, tf1, o1_bv);
    computeBV(s2, tf2, o2_bv);
    return !BoxPairTranslationalCCD<S>::IsDisjoint(
        o1_bv, s1_displacement, o2_bv, toc, request.zero_movement_tolerance);
  }
};

template <typename S_>
template <typename Shape1, typename Shape2>
void ShapePairTranslationalCollisionSolver<S_>::RunShapePair(
    const Shape1& s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement, const Shape2& s2,
    const Transform3<S>& tf2, const ContinuousContactMeta<S>& contact_meta,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Enough collision result, do NOT continue
  if (result.TerminationConditionSatisfied(request.num_max_contacts)) {
    return;
  }

  // Request type
  const TimeOfCollisionRequestType request_type =
      (request.request_type == TimeOfCollisionRequestType::kBoxApproximate &&
       contact_meta.external_toc_valid())
          ? TimeOfCollisionRequestType::kNotRequested
          : request.request_type;

  // Invoke interface
  Interval<S> toc;
  bool is_intersect =
      ShapePairTranslationalCollisionImpl<S, Shape1, Shape2>::RunIntersect(
          s1, tf1, s1_displacement, s2, tf2, request, request_type, toc);
  if (!is_intersect) {
    return;
  }

  // Push the contact
  ContinuousCollisionContact<S> contact;
  contact_meta.writeToContact(contact, toc);
  result.AddContact(std::move(contact));
}

template <typename S_>
template <typename Shape1, typename Shape2>
void ShapePairTranslationalCollisionSolver<S_>::RunShapePair(
    const Shape1* s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement, const Shape2* s2,
    const Transform3<S>& tf2, const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // Write the meta just as shape pair
  ContinuousContactMeta<S> contact_meta;
  contact_meta.o1 = s1;
  contact_meta.o2 = s2;
  contact_meta.b1 = ContinuousContactMeta<S>::kNone;
  contact_meta.b2 = ContinuousContactMeta<S>::kNone;
  contact_meta.external_box_toc.lower_bound = -1;
  contact_meta.external_box_toc.upper_bound = -1;

  // Invoke interface
  RunShapePair(*s1, tf1, s1_displacement, *s2, tf2, contact_meta, request,
               result);
}

template <typename S_>
template <typename Shape>
void ShapePairTranslationalCollisionSolver<S_>::RunShapeSimplex(
    const Shape& s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement, const Simplex<S>& s2,
    const Transform3<S>& tf2, const ContinuousContactMeta<S>& contact_meta,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  if (s2.is_triangle()) {
    const TriangleP<S> triangle(s2[0], s2[1], s2[2]);
    RunShapePair<Shape, TriangleP<S>>(s1, tf1, s1_displacement, triangle, tf2,
                                      contact_meta, request, result);
  } else {
    assert(s2.is_tetrahedron());
    const Tetrahedron<S> tetrahedron(s2[0], s2[1], s2[2], s2[3]);
    RunShapePair<Shape, Tetrahedron<S>>(s1, tf1, s1_displacement, tetrahedron,
                                        tf2, contact_meta, request, result);
  }
}

template <typename S_>
void ShapePairTranslationalCollisionSolver<S_>::RunSimplexPair(
    const Simplex<S>& s1, const Transform3<S>& tf1,
    const TranslationalDisplacement<S>& s1_displacement, const Simplex<S>& s2,
    const Transform3<S>& tf2, const ContinuousContactMeta<S>& contact_meta,
    const ContinuousCollisionRequest<S>& request,
    ContinuousCollisionResult<S>& result) {
  // No primitive
  if (s1.is_triangle() && s2.is_triangle()) {
    const TriangleP<S> tri1(s1[0], s1[1], s1[2]);
    const TriangleP<S> tri2(s2[0], s2[1], s2[2]);
    RunShapePair<TriangleP<S>, TriangleP<S>>(
        tri1, tf1, s1_displacement, tri2, tf2, contact_meta, request, result);
  } else if (s1.is_triangle() && s2.is_tetrahedron()) {
    const TriangleP<S> tri(s1[0], s1[1], s1[2]);
    const Tetrahedron<S> tet(s2[0], s2[1], s2[2], s2[3]);
    RunShapePair<TriangleP<S>, Tetrahedron<S>>(
        tri, tf1, s1_displacement, tet, tf2, contact_meta, request, result);
  } else if (s1.is_tetrahedron() && s2.is_triangle()) {
    const Tetrahedron<S> tet(s1[0], s1[1], s1[2], s1[3]);
    const TriangleP<S> tri(s2[0], s2[1], s2[2]);
    RunShapePair<Tetrahedron<S>, TriangleP<S>>(
        tet, tf1, s1_displacement, tri, tf2, contact_meta, request, result);
  } else {
    assert(s1.is_tetrahedron() && s2.is_tetrahedron());
    const Tetrahedron<S> tet1(s1[0], s1[1], s1[2], s1[3]);
    const Tetrahedron<S> tet2(s2[0], s2[1], s2[2], s2[3]);
    RunShapePair<Tetrahedron<S>, Tetrahedron<S>>(
        tet1, tf1, s1_displacement, tet2, tf2, contact_meta, request, result);
  }
}

}  // namespace detail
}  // namespace fcl