#pragma once

#include "fcl/narrowphase/detail/traversal/collision/intersect.h"

namespace fcl {
namespace detail {

template <typename S>
void ContactMeta<S>::reset() {
  o1 = nullptr;
  o2 = nullptr;
  b1 = kNone;
  b2 = kNone;
  reverse_normal = false;
  reverse_o1_and_o2 = false;
}

template <typename S>
void ContactMeta<S>::writeToContact(const ContactPoint<S>& contact_point,
                                    Contact<S>& contact) const {
  writeToContact(contact);
  contact.pos = contact_point.pos;
  contact.penetration_depth = contact_point.penetration_depth;
  contact.normal = contact_point.normal;
  if (reverse_normal) contact.normal *= -1;
}

template <typename S>
void ContactMeta<S>::writeToContact(Contact<S>& contact) const {
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
}

template <typename S_>
template <typename Shape1, typename Shape2>
void ShapePairIntersectSolver<S_>::ShapeIntersect(
    const Shape1& s1, const Transform3<S_>& tf1, const Shape2& s2,
    const Transform3<S_>& tf2, const CollisionRequest<S_>& request,
    const ContactMeta<S>& contact_meta, CollisionResult<S_>& result) const {
  static_assert(std::is_same<S_, typename Shape1::S>::value,
                "scalar type must match");
  static_assert(std::is_same<S_, typename Shape2::S>::value,
                "scalar type must match");
  // Enough collision result, do NOT continue
  if (result.numContacts() >= request.maxNumContacts()) {
    return;
  }

  // Not occupied
  if ((!s1.isOccupied()) || (!s2.isOccupied())) {
    return;
  }

  // Do NOT need penetration
  if (!request.isPenetrationEnabled()) {
    const bool is_intersect =
        gjk_solver->shapeIntersect(s1, tf1, s2, tf2, nullptr);
    if (is_intersect) {
      assert(result.numContacts() < request.maxNumContacts());
      Contact<S> contact;
      contact_meta.writeToContact(contact);
      result.addContact(contact);
    }

    // Done with no penetration case
    return;
  }

  assert(request.isPenetrationEnabled());
  CollisionPenetrationContactData<S> contacts;
  const bool is_intersect =
      gjk_solver->shapeIntersect(s1, tf1, s2, tf2, &contacts);
  if (!is_intersect) {
    // No collision, direct return
    return;
  }

  assert(result.numContacts() < request.maxNumContacts());
  const std::size_t free_space =
      request.maxNumContacts() - result.numContacts();
  std::size_t num_adding_contacts{0};

  // If the free space is not enough to add all the new contacts, we add
  // contacts in descent order of penetration depth.
  if (free_space < contacts.size()) {
    std::partial_sort(contacts.begin(), contacts.begin() + free_space,
                      contacts.end(),
                      std::bind(comparePenDepth<S>, std::placeholders::_2,
                                std::placeholders::_1));
    num_adding_contacts = free_space;
  } else {
    num_adding_contacts = contacts.size();
  }

  // Add the contacts
  for (std::size_t i = 0; i < num_adding_contacts; ++i) {
    Contact<S> contact_i;
    contact_meta.writeToContact(contacts[i], contact_i);
    result.addContact(contact_i);
  }
}

template <typename S_>
template <typename Shape1, typename Shape2>
void ShapePairIntersectSolver<S_>::ShapeIntersect(
    const Shape1* s1, const Transform3<S_>& tf1, const Shape2* s2,
    const Transform3<S_>& tf2, const CollisionRequest<S_>& request,
    CollisionResult<S_>& result) const {
  // Write the meta just as shape pair
  ContactMeta<S> contact_meta;
  contact_meta.o1 = s1;
  contact_meta.o2 = s2;
  contact_meta.b1 = ContactMeta<S>::kNone;
  contact_meta.b2 = ContactMeta<S>::kNone;
  ShapeIntersect(*s1, tf1, *s2, tf2, request, contact_meta, result);
}

template <typename S_>
template <typename Shape>
void ShapePairIntersectSolver<S_>::ShapeSimplexIntersect(
    const Shape& s1, const Transform3<S_>& tf1, const Simplex<S_>& s2,
    const Transform3<S_>& tf2, const CollisionRequest<S_>& request,
    const ContactMeta<S>& contact_meta, CollisionResult<S_>& result) const {
  static_assert(std::is_same<S_, typename Shape::S>::value,
                "scalar type must match");
  // Enough collision result, do NOT continue
  if (result.numContacts() >= request.maxNumContacts()) {
    return;
  }

  // Not occupied, only for shape
  if (!s1.isOccupied()) {
    return;
  }

  if (s2.is_triangle()) {
    if (!request.isPenetrationEnabled()) {
      bool intersect = gjk_solver->shapeTriangleIntersect(s1, tf1, s2[0], s2[1],
                                                          s2[2], tf2, nullptr);
      if (intersect) {
        assert(result.numContacts() < request.maxNumContacts());
        Contact<S> this_contact;
        contact_meta.writeToContact(this_contact);
        result.addContact(std::move(this_contact));
      }
    } else {
      // The contact result
      CollisionPenetrationContactData<S> penetration_data;
      bool intersect = gjk_solver->shapeTriangleIntersect(
          s1, tf1, s2[0], s2[1], s2[2], tf2, &penetration_data);
      if (intersect) {
        assert(result.numContacts() < request.maxNumContacts());
        assert(penetration_data.size() == 1u);
        const auto& this_penetration = penetration_data[0];
        Contact<S> this_contact;
        contact_meta.writeToContact(this_penetration, this_contact);
        result.addContact(std::move(this_contact));
      }
    }
  } else {
    if (!request.isPenetrationEnabled()) {
      bool intersect = gjk_solver->shapeTetrahedronIntersect(
          s1, tf1, s2[0], s2[1], s2[2], s2[3], tf2, nullptr);
      if (intersect) {
        assert(result.numContacts() < request.maxNumContacts());
        Contact<S> this_contact;
        contact_meta.writeToContact(this_contact);
        result.addContact(std::move(this_contact));
      }
    } else {
      // The contact result
      CollisionPenetrationContactData<S> penetration_data;
      bool intersect = gjk_solver->shapeTetrahedronIntersect(
          s1, tf1, s2[0], s2[1], s2[2], s2[3], tf2, &penetration_data);
      if (intersect) {
        assert(result.numContacts() < request.maxNumContacts());
        assert(penetration_data.size() == 1u);
        const auto& this_penetration = penetration_data[0];
        Contact<S> this_contact;
        contact_meta.writeToContact(this_penetration, this_contact);
        result.addContact(std::move(this_contact));
      }
    }
  }
}

template <typename S_>
void ShapePairIntersectSolver<S_>::trianglePairIntersect(
    const Simplex<S_>& s1, const Transform3<S_>& tf1, const Simplex<S_>& s2,
    const Transform3<S_>& tf2, const Matrix3<S>& rotation_2to1,
    const Vector3<S>& translation_2in1, const CollisionRequest<S_>& request,
    const ContactMeta<S>& contact_meta, CollisionResult<S_>& result) const {
  FCL_UNUSED(tf2);
  assert(s1.is_triangle() && s2.is_triangle());
  if (!request.isPenetrationEnabled()) {
    // only interested in collision or not
    if (Intersect<S>::intersect_Triangle(s1[0], s1[1], s1[2], s2[0], s2[1],
                                         s2[2], rotation_2to1,
                                         translation_2in1)) {
      assert(result.numContacts() < request.maxNumContacts());
      Contact<S> this_contact;
      contact_meta.writeToContact(this_contact);
      result.addContact(this_contact);
    }

    return;
  }

  // need compute the contact information
  assert(request.isPenetrationEnabled());
  S penetration;
  Vector3<S> normal;
  unsigned int n_contacts;
  Vector3<S> contacts[2];

  if (Intersect<S>::intersect_Triangle(
          s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], rotation_2to1,
          translation_2in1, contacts, &n_contacts, &penetration, &normal)) {
    if (request.maxNumContacts() < result.numContacts() + n_contacts) {
      n_contacts = (request.maxNumContacts() > result.numContacts())
                       ? (request.maxNumContacts() - result.numContacts())
                       : 0;
    }

    for (unsigned int i = 0; i < n_contacts; ++i) {
      Contact<S> contact_i;
      contact_meta.writeToContact(contact_i);
      const bool reverse_normal_i = contact_meta.reverse_normal;

      // Into world frame
      contact_i.pos = tf1 * contacts[i];
      contact_i.penetration_depth = penetration;
      contact_i.normal = tf1.linear() * normal;
      if (reverse_normal_i) contact_i.normal *= -1;
      result.addContact(contact_i);
    }
  }
}

template <typename S_>
void ShapePairIntersectSolver<S_>::SimplexIntersect(
    const Simplex<S_>& s1, const Transform3<S_>& tf1, const Simplex<S_>& s2,
    const Transform3<S_>& tf2, const Matrix3<S>& rotation_2to1,
    const Vector3<S>& translation_2in1, const CollisionRequest<S_>& request,
    const ContactMeta<S>& contact_meta, CollisionResult<S_>& result) const {
  // Enough collision result, do NOT continue
  if (result.numContacts() >= request.maxNumContacts()) {
    return;
  }

  // Only triangle-triangle is primitive
  if (s1.is_triangle() && s2.is_triangle()) {
    trianglePairIntersect(s1, tf2, s2, tf2, rotation_2to1, translation_2in1,
                          request, contact_meta, result);
  } else if (s1.is_triangle() && s2.is_tetrahedron()) {
    const TriangleP<S> tri(s1[0], s1[1], s1[2]);
    const Tetrahedron<S> tet(s2[0], s2[1], s2[2], s2[3]);
    ShapeIntersect<TriangleP<S>, Tetrahedron<S>>(tri, tf1, tet, tf2, request,
                                                 contact_meta, result);
  } else if (s1.is_tetrahedron() && s2.is_triangle()) {
    const Tetrahedron<S> tet(s1[0], s1[1], s1[2], s1[3]);
    const TriangleP<S> tri(s2[0], s2[1], s2[2]);
    ShapeIntersect<Tetrahedron<S>, TriangleP<S>>(tet, tf1, tri, tf2, request,
                                                 contact_meta, result);
  } else {
    assert(s1.is_tetrahedron() && s2.is_tetrahedron());
    const Tetrahedron<S> tet1(s1[0], s1[1], s1[2], s1[3]);
    const Tetrahedron<S> tet2(s2[0], s2[1], s2[2], s2[3]);
    ShapeIntersect<Tetrahedron<S>, Tetrahedron<S>>(
        tet1, tf1, tet2, tf2, request, contact_meta, result);
  }
}

}  // namespace detail
}  // namespace fcl