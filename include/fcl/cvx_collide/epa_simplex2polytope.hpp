//
// Created by wei on 22-6-5.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
void extractTouchingPoint(const MinkowskiDiff<T>& shape,
                          const Vector3<T>& direction,
                          Vector3<T>* p0_if_touching,
                          Vector3<T>* p1_if_touching) {
  // In touching, two objects should have zero separation/penetration distance,
  // and the contact should be a point contact. However, in some situation the
  // contact can be a line/face contact EVEN AT TOUCHING, for example two boxes
  // touch each other via faces.
  // Previously, we compute the point on two objects as:
  //     p0 = shape.support0(direction);
  //     p1 = shape.support1(-direction);
  // However, in the box example above p0 and p1 computed above might
  // NOT BE SAME (ROUGHLY). Actually, their distance can be arbitrarily large if
  // the direction is perpendicular to the touching face, as any point on the
  // touching faces is a valid output of the support function.
  //
  // Solution:
  // As p0 and p1 are both expressed in object0 frame and the should be the
  // same, we may only p0 and set p1 = p0, or use the middle point
  if ((p0_if_touching == nullptr) && (p1_if_touching == nullptr)) {
    return;
  }

  // Compute the middle point
  const Vector3<T> p_middle =
      (shape.support0(direction) + shape.support1(-direction)) / 2;
  if (p0_if_touching) *p0_if_touching = p_middle;
  if (p1_if_touching) *p1_if_touching = p_middle;

  // old impl
  // if (p0_if_touching) *p0_if_touching = shape.support0(direction);
  // if (p1_if_touching) *p1_if_touching = shape.support1(-direction);
}

template <typename T>
Simplex2PolytopeStatus simplexToPolytope(const GJKSimplex<T>& simplex,
                                         const MinkowskiDiff<T>& shape,
                                         Polytope<T>& polytope,
                                         Vector3<T>* p0_if_touching,
                                         Vector3<T>* p1_if_touching,
                                         T touching_threshold) {
  // Check non-of them are close to origin, else just report touching
  const auto simplex_size = simplex.rank;
  const T touching_threshold_square = touching_threshold * touching_threshold;
  for (int i = 0; i < simplex_size; i++) {
    const auto& v_i = simplex.vertices[i];
    const bool point_touching_i =
        v_i.vertex.squaredNorm() <= touching_threshold_square;
    if (point_touching_i) {
      extractTouchingPoint(shape, v_i.direction, p0_if_touching,
                           p1_if_touching);
      return Simplex2PolytopeStatus::Touching;
    }
  }

  // No-point touching case
  if (simplex_size == 4) {
    const auto& a = simplex.vertices[0];
    const auto& b = simplex.vertices[1];
    const auto& c = simplex.vertices[2];
    const auto& d = simplex.vertices[3];
    return simplexToPolytope4(a, b, c, d, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  } else if (simplex_size == 3) {
    const auto& a = simplex.vertices[0];
    const auto& b = simplex.vertices[1];
    const auto& c = simplex.vertices[2];
    return simplexToPolytope3(a, b, c, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  } else if (simplex_size == 2) {
    const auto& a = simplex.vertices[0];
    const auto& b = simplex.vertices[1];
    return simplexToPolytope2(a, b, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  } else {
    // Only one point?
    const auto& a = simplex.vertices[0];
    const Vector3<T>& direction = a.direction;
    if (p0_if_touching) *p0_if_touching = shape.support0(direction);
    if (p1_if_touching) *p1_if_touching = shape.support1(-direction);
    return Simplex2PolytopeStatus::Touching;
  }
}

template <typename T>
Simplex2PolytopeStatus simplexToPolytope4(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiffVertex<T>& d,
    const MinkowskiDiff<T>& shape, Polytope<T>& polytope,
    Vector3<T>* p0_if_touching, Vector3<T>* p1_if_touching,
    T touching_threshold) {
  // The origin
  const Vector3<T> o = Vector3<T>::Zero();

  // Goes to abc
  if (pointToPlaneDistance(o, a.vertex, b.vertex, c.vertex) <
      touching_threshold) {
    return simplexToPolytope3(a, b, c, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  }

  if (pointToPlaneDistance(o, a.vertex, c.vertex, d.vertex) <
      touching_threshold) {
    return simplexToPolytope3(a, c, d, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  }

  if (pointToPlaneDistance(o, a.vertex, b.vertex, d.vertex) <
      touching_threshold) {
    return simplexToPolytope3(a, b, d, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  }

  if (pointToPlaneDistance(o, b.vertex, c.vertex, d.vertex) <
      touching_threshold) {
    return simplexToPolytope3(b, c, d, shape, polytope, p0_if_touching,
                              p1_if_touching, touching_threshold);
  }

  // No touching, simply create tetrahedron
  return formNewTetrahedronPolytope(a, b, c, d, polytope);
}

template <typename T>
Simplex2PolytopeStatus simplexToPolytope3(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiff<T>& shape,
    Polytope<T>& polytope, Vector3<T>* p0_if_touching,
    Vector3<T>* p1_if_touching, T touching_threshold) {
  const Vector3<T> ab = b.vertex - a.vertex;
  const Vector3<T> ac = c.vertex - a.vertex;
  Vector3<T> abc_normal = ab.cross(ac);

  // Make a unit vector
  if (abc_normal.squaredNorm() <= 0.0) return Simplex2PolytopeStatus::Failed;
  abc_normal.normalize();

  // On normal/-normal
  const Vector3<T>& direction_0 = abc_normal;
  const Vector3<T> d0 = shape.support(direction_0);
  const T d0_to_plane = pointToPlaneDistance(d0, a.vertex, b.vertex, c.vertex);
  const Vector3<T>& direction_1 = -abc_normal;
  const Vector3<T> d1 = shape.support(direction_1);
  const T d1_to_plane = pointToPlaneDistance(d1, a.vertex, b.vertex, c.vertex);
  if (d0_to_plane <= touching_threshold) {
    extractTouchingPoint(shape, direction_0, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  } else if (d1_to_plane <= touching_threshold) {
    extractTouchingPoint(shape, direction_1, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  }

  if (d0_to_plane > d1_to_plane) {
    MinkowskiDiffVertex<T> d0_vertex;
    d0_vertex.vertex = d0;
    d0_vertex.direction = direction_0;
    return formNewTetrahedronPolytope(a, b, c, d0_vertex, polytope);
  } else {
    MinkowskiDiffVertex<T> d1_vertex;
    d1_vertex.vertex = d1;
    d1_vertex.direction = direction_1;
    return formNewTetrahedronPolytope(a, b, c, d1_vertex, polytope);
  }
}

template <typename T>
Simplex2PolytopeStatus formNewTetrahedronPolytope(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiffVertex<T>& d,
    Polytope<T>& polytope) {
  // Re-init the polytope
  auto face_capacity = polytope.face_capacity();
  polytope.Reset(face_capacity);

  // 4 vertices of the tetrahedron
  std::array<PolytopeVertex<T>*, 4> v;
  v[0] = polytope.AddNewVertex(a);
  v[1] = polytope.AddNewVertex(b);
  v[2] = polytope.AddNewVertex(c);
  v[3] = polytope.AddNewVertex(d);

  // 6 edges
  std::array<PolytopeEdge<T>*, 6> e;
  e[0] = polytope.AddNewEdge(v[0], v[1]);
  e[1] = polytope.AddNewEdge(v[1], v[2]);
  e[2] = polytope.AddNewEdge(v[2], v[0]);
  e[3] = polytope.AddNewEdge(v[3], v[0]);
  e[4] = polytope.AddNewEdge(v[3], v[1]);
  e[5] = polytope.AddNewEdge(v[3], v[2]);

  // 4 faces
  std::array<PolytopeFace<T>*, 4> f;
  f[0] = polytope.AddNewFace(e[0], e[1], e[2]);
  f[1] = polytope.AddNewFace(e[3], e[4], e[0]);
  f[2] = polytope.AddNewFace(e[4], e[5], e[1]);
  f[3] = polytope.AddNewFace(e[5], e[3], e[2]);

  // Only check at the face level
  for (auto& face_i : f) {
    if (face_i == nullptr) return Simplex2PolytopeStatus::Failed;
  }

  // OK
  return Simplex2PolytopeStatus::OK;
}

template <typename T>
Simplex2PolytopeStatus simplexToPolytope2(const MinkowskiDiffVertex<T>& a,
                                          const MinkowskiDiffVertex<T>& b,
                                          const MinkowskiDiff<T>& shape,
                                          Polytope<T>& polytope,
                                          Vector3<T>* p0_if_touching,
                                          Vector3<T>* p1_if_touching,
                                          T touching_threshold) {
  // Step 1: search for a point on the shape that is outside a/b
  //         if not exists, then report touching
  const T distance_threshold = touching_threshold;
  const T distance_threshold_squared = touching_threshold * touching_threshold;
  const Vector3<T> a_to_b = (b.vertex - a.vertex);
  if (a_to_b.squaredNorm() < distance_threshold_squared)
    return Simplex2PolytopeStatus::Failed;
  const Vector3<T> a_to_b_unit = a_to_b.normalized();

  // v0 is any vector that perpendicular to a_to_b_unit
  Vector3<T> d_init(a_to_b_unit[1], -a_to_b_unit[0], 0);
  bool d_init_valid = (d_init.squaredNorm() > distance_threshold_squared);

  // Search for other options in case v0 is zero
  {
    // Attempt 1
    if (!d_init_valid) {
      d_init = Vector3<T>(a_to_b_unit[2], 0, -a_to_b_unit[0]);
      d_init_valid = (d_init.squaredNorm() > distance_threshold_squared);
    }

    // Attempt 2
    if (!d_init_valid) {
      d_init = Vector3<T>(0, a_to_b_unit[2], -a_to_b_unit[1]);
      d_init_valid = (d_init.squaredNorm() > distance_threshold_squared);
    }
  }

  // Should be valid
  assert(d_init_valid);
  if (!d_init_valid) return Simplex2PolytopeStatus::Failed;
  d_init.normalize();

  // Start searching
  constexpr int div_2pi_by_n = 72;
  Vector3<T> d_current = d_init;
  Vector3<T> d_not_on_ab;
  Vector3<T> v_not_on_ab;
  bool v_not_on_ab_valid = false;
  for (auto i = 0; i < div_2pi_by_n; i++) {
    // Is v current not on ab?
    const Vector3<T> v_current = shape.support(d_current);
    const T distance = pointToLineDistance(v_current, a.vertex, b.vertex);
    if (distance > distance_threshold) {
      v_not_on_ab = v_current;
      d_not_on_ab = d_current;
      v_not_on_ab_valid = true;
      break;
    }

    // The next one
    constexpr T pi_value = 3.145926;
    constexpr T delta_vector_length = T(2.0) * pi_value / div_2pi_by_n;
    const Vector3<T> delta_vector = a_to_b_unit.cross(v_current);
    d_current = d_current + delta_vector_length * delta_vector;
    d_current.normalize();
  }

  // Everything on ab, touching
  if (!v_not_on_ab_valid) {
    extractTouchingPoint(shape, d_init, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  }

  // Now we have a point not on ab
  MinkowskiDiffVertex<T> v0;
  v0.vertex = v_not_on_ab;
  v0.direction = d_not_on_ab;

  // Get second support point in opposite direction than d0
  MinkowskiDiffVertex<T> v1;
  v1.direction = -v0.direction;
  v1.vertex = shape.support(v1.direction);
  if (pointToLineDistance(v1.vertex, a.vertex, b.vertex) < distance_threshold) {
    extractTouchingPoint(shape, v1.direction, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  }

  // Next will be in direction of normal of triangle a, v0, v1
  const Vector3<T> a_to_v0 = v0.vertex - a.vertex;
  const Vector3<T> a_to_v1 = v1.vertex - a.vertex;
  const Vector3<T> a_v0_v1_normal = a_to_v0.cross(a_to_v1);
  if (a_v0_v1_normal.squaredNorm() < distance_threshold_squared)
    return Simplex2PolytopeStatus::Failed;

  // v2
  MinkowskiDiffVertex<T> v2;
  v2.direction = a_v0_v1_normal.normalized();
  v2.vertex = shape.support(v2.direction);
  if (pointToLineDistance(v2.vertex, a.vertex, b.vertex) < distance_threshold) {
    extractTouchingPoint(shape, v2.direction, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  }

  // Last point
  MinkowskiDiffVertex<T> v3;
  v3.direction = -v2.direction;
  v3.vertex = shape.support(v3.direction);
  if (pointToLineDistance(v3.vertex, a.vertex, b.vertex) < distance_threshold) {
    extractTouchingPoint(shape, v3.direction, p0_if_touching, p1_if_touching);
    return Simplex2PolytopeStatus::Touching;
  }

  // No touching, we have constructed enclosing polytopes
  {
    // Re-init the polytope
    auto face_capacity = polytope.face_capacity();
    polytope.Reset(face_capacity);

    // 6 vertices of two tetrahedron
    std::array<PolytopeVertex<T>*, 6> v;
    v[0] = polytope.AddNewVertex(a);
    v[1] = polytope.AddNewVertex(v0);
    v[2] = polytope.AddNewVertex(b);
    v[3] = polytope.AddNewVertex(v1);
    v[4] = polytope.AddNewVertex(v2);
    v[5] = polytope.AddNewVertex(v3);

    // 12 edges
    std::array<PolytopeEdge<T>*, 12> e;
    e[0] = polytope.AddNewEdge(v[0], v[1]);
    e[1] = polytope.AddNewEdge(v[1], v[2]);
    e[2] = polytope.AddNewEdge(v[2], v[3]);
    e[3] = polytope.AddNewEdge(v[3], v[0]);

    e[4] = polytope.AddNewEdge(v[4], v[0]);
    e[5] = polytope.AddNewEdge(v[4], v[1]);
    e[6] = polytope.AddNewEdge(v[4], v[2]);
    e[7] = polytope.AddNewEdge(v[4], v[3]);

    e[8] = polytope.AddNewEdge(v[5], v[0]);
    e[9] = polytope.AddNewEdge(v[5], v[1]);
    e[10] = polytope.AddNewEdge(v[5], v[2]);
    e[11] = polytope.AddNewEdge(v[5], v[3]);

    // 8 faces
    std::array<PolytopeFace<T>*, 8> f;
    f[0] = polytope.AddNewFace(e[4], e[5], e[0]);
    f[1] = polytope.AddNewFace(e[5], e[6], e[1]);
    f[2] = polytope.AddNewFace(e[6], e[7], e[2]);
    f[3] = polytope.AddNewFace(e[7], e[4], e[3]);

    f[4] = polytope.AddNewFace(e[8], e[9], e[0]);
    f[5] = polytope.AddNewFace(e[9], e[10], e[1]);
    f[6] = polytope.AddNewFace(e[10], e[11], e[2]);
    f[7] = polytope.AddNewFace(e[11], e[8], e[3]);

    // Only check at the face level
    for (auto& face_i : f) {
      if (face_i == nullptr) return Simplex2PolytopeStatus::Failed;
    }
  }

  // Done
  return Simplex2PolytopeStatus::OK;
}

}  // namespace cvx_collide
}  // namespace fcl