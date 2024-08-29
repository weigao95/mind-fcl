//
// Created by wei on 22-6-18.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
GJK_Status GJK<T>::Evaluate(
    const MinkowskiDiff<T>& shape, GJKSimplex<T>& simplex,
    const Vector3<T>& guess,
    MinSeparationDistanceOutput* min_distance_output_if_separated) const {
  // The initial guess
  Vector3<T> direction = guess;
  if (direction.squaredNorm() <= 0.0) direction = Vector3<T>::UnitX();
  direction.normalize();

  // Handle of separation
  auto process_separated_vertex =
      [this, &shape, &min_distance_output_if_separated,
       &simplex](const MinkowskiDiffVertex<T>& separated_v) -> GJK_Status {
    // Do not require distance, just return
    if (min_distance_output_if_separated == nullptr)
      return GJK_Status::Separated;

    // Assign output
    min_distance_output_if_separated->is_certification_vertex_valid = true;
    min_distance_output_if_separated->vertex_certify_separation = separated_v;

    // Init the simplex as the separation point
    simplex.reset();
    simplex.AddVertex(separated_v);

    // Compute the output
    std::pair<Vector3<T>, Vector3<T>> min_distance_raw_output;
    auto distance_ok = findMinimumDistancePointsWithSeparatedVertexInit(
        shape, simplex, min_distance_raw_output);

    // Assign
    min_distance_output_if_separated->is_separation_point_valid = distance_ok;
    if (distance_ok) {
      min_distance_output_if_separated->p0_if_separated =
          min_distance_raw_output.first;
      min_distance_output_if_separated->p1_if_separated =
          min_distance_raw_output.second;
    }

    return GJK_Status::Separated;
  };

  auto process_undetermined_stop =
      [this, &shape, &min_distance_output_if_separated,
       &simplex](GJK_Status return_status) -> GJK_Status {
    if (min_distance_output_if_separated != nullptr) {
      // No separation certificate point
      min_distance_output_if_separated->is_certification_vertex_valid = false;

      // Compute the distance
      std::pair<Vector3<T>, Vector3<T>> min_distance_raw_output;
      auto distance_ok = extractSeparationPointTrySubSimplex(
          shape, simplex, min_distance_raw_output);
      min_distance_output_if_separated->is_separation_point_valid = distance_ok;
      if (distance_ok) {
        min_distance_output_if_separated->p0_if_separated =
            min_distance_raw_output.first;
        min_distance_output_if_separated->p1_if_separated =
            min_distance_raw_output.second;
      }
    }
    return return_status;
  };

  // Init the simplex first vertex
  MinkowskiDiffVertex<T> vertex = shape.supportVertex(direction);
  simplex.reset();
  simplex.AddVertex(vertex);
  assert(simplex.n_vertices() == 1);

  // The vertex is almost origin
  const T tolerance_squared = tolerance_ * tolerance_;
  if (vertex.vertex.squaredNorm() <= tolerance_squared) {
    // Just use this vertex, but return
    return GJK_Status::Intersect;
  } else if (vertex.vertex.dot(direction) < 0) {
    return process_separated_vertex(vertex);
  }

  // Revert the direction and start the iteration
  direction *= -1;
  std::size_t iteration = 0;
  while (iteration < max_iterations_) {
    // Update the index
    iteration += 1;

    // Obtain the support point
    vertex = shape.supportVertex(direction);

    // Check if farthest point in Minkowski difference in direction dir
    // isn't somewhere before origin (the test on negative dot product)
    // - because if it is, objects are not intersecting at all.
    if (vertex.vertex.dot(direction) < 0) {
      return process_separated_vertex(vertex);
    }

    // Should not be the same as old ones
    for (auto j = 0; j < simplex.n_vertices(); j++) {
      const Vector3<T>& vertex_j = simplex.vertices[j].vertex;
      if ((vertex_j - vertex.vertex).squaredNorm() < tolerance_squared) {
        return process_undetermined_stop(GJK_Status::ConvergeNoProgress);
      }
    }

    // The new vertex is almost origin
    if (vertex.vertex.squaredNorm() <= tolerance_squared) {
      // Just use this vertex, but return
      simplex.AddVertex(vertex);
      return GJK_Status::Intersect;
    }

    // Add the simplex to vertex
    simplex.AddVertex(vertex);

    // Do projection
    ProjectionStatus status = simplexProjection(simplex, direction);
    if (status == ProjectionStatus::Failed) {
      return GJK_Status::Failed;
    } else if (status == ProjectionStatus::Intersect) {
      return GJK_Status::Intersect;
    } else if (status == ProjectionStatus::FailedTetrahedronZeroVolume) {
      // Tetrahedron are obtained by extending a triangle on its normal
      // The triangle must have positive area, thus the new vertex must
      // lie on the plane of old triangle, and we converge
      return process_undetermined_stop(GJK_Status::ConvergeNoProgress);
    }
  }

  return process_undetermined_stop(GJK_Status::IterationLimit);
}

template <typename T>
typename GJK<T>::ProjectionStatus GJK<T>::simplexProjection(
    GJKSimplex<T>& simplex, Vector3<T>& direction) const {
  // assert: no duplicate vertex (|a - b| > tolerance)
  //         no vertex closet to origin (|a| > tolerance)
  assert(simplex.n_vertices() >= 2);
  if (simplex.n_vertices() == 2) {
    return simplexProjection2(simplex, direction);
  }

  if (simplex.n_vertices() == 3) {
    return simplexProjection3(simplex, direction);
  }

  // Should be 4-vertex
  assert(simplex.n_vertices() == 4);
  return simplexProjection4(simplex, direction);
}

template <typename T>
typename GJK<T>::ProjectionStatus GJK<T>::simplexProjection2(
    GJKSimplex<T>& simplex, Vector3<T>& direction) const {
  // A is the newly added vertex
  MinkowskiDiffVertex<T> B = simplex.vertices[0];
  MinkowskiDiffVertex<T> A = simplex.vertices[1];

  // Now we can use reference
  const Vector3<T>& a = A.vertex;
  const Vector3<T>& b = B.vertex;
  const Vector3<T> a_to_b = b - a;
  const Vector3<T>& o_to_a = a;
  const T ao_dot_ab = -a_to_b.dot(o_to_a);

  // Check if the origin is on AB segment
  // assert o_to_b.norm() >= tolerance_
  // assert o_to_a.norm() >= tolerance_
  // This is oa_cross_ab, which equals ab_cross_ao
  const Vector3<T> ab_cross_ao = o_to_a.cross(a_to_b);

  // Here we use a very strict 0.0 tolerance on segment
  if (ao_dot_ab > 0 && ab_cross_ao.squaredNorm() <= 0.0) {
    return ProjectionStatus::Intersect;
  }

  // Not on AB segment
  if (ao_dot_ab <= 0) {
    // The closest point is A
    simplex.reset();
    simplex.AddVertex(A);
    direction = -o_to_a.normalized();
    return ProjectionStatus::Continue;
  } else {
    // The closest point is on AB segment
    // Do not touch the AB simplex
    direction = ab_cross_ao.cross(a_to_b);
    assert(o_to_a.dot(direction) <= 0);
    direction.normalize();
    return ProjectionStatus::Continue;
  }
}

template <typename T>
typename GJK<T>::ProjectionStatus GJK<T>::simplexProjection3(
    GJKSimplex<T>& simplex, Vector3<T>& direction) const {
  // A is the newly added vertex
  MinkowskiDiffVertex<T> C = simplex.vertices[0];
  MinkowskiDiffVertex<T> B = simplex.vertices[1];
  MinkowskiDiffVertex<T> A = simplex.vertices[2];

  const Vector3<T>& a = A.vertex;
  const Vector3<T>& b = B.vertex;
  const Vector3<T>& c = C.vertex;

  // The case of a
  const Vector3<T> a_to_b = b - a;
  const Vector3<T> a_to_c = c - a;
  const Vector3<T>& o_to_a = a;
  const bool a_separate_b_o = o_to_a.dot(a_to_b) >= 0;
  const bool a_separate_c_o = o_to_a.dot(a_to_c) >= 0;
  if (a_separate_b_o && a_separate_c_o) {
    simplex.reset();
    simplex.AddVertex(A);
    direction = -o_to_a;
    direction.normalize();
    return ProjectionStatus::Continue;
  }

  // We know that |a_to_b| > tolerance, |a_to_c| > tolerance
  // Thus, the normal is not a zero vector
  const Vector3<T> ab_cross_ac = a_to_b.cross(a_to_c);
  const Vector3<T>& abc_normal = ab_cross_ac;

  // This guarantee a_to_b.dot(ac_normal_in_abc) <= 0
  const Vector3<T> ac_normal_in_abc = abc_normal.cross(a_to_c);
  assert(a_to_b.dot(ac_normal_in_abc) <= 0);
  if (o_to_a.dot(ac_normal_in_abc) <= 0) {
    // ON AC, simplex.vertices[0] == C
    simplex.vertices[1] = A;
    simplex.rank = 2;

    // direction is the normal of AC on OAC plane
    const Vector3<T> ac_cross_ao = o_to_a.cross(a_to_c);
    direction = ac_cross_ao.cross(a_to_c);
    direction.normalize();
    assert(o_to_a.dot(direction) <= 0);
    return ProjectionStatus::Continue;
  }

  // This guarantee a_to_c.dot(ab_normal_in_abc) <= 0
  const Vector3<T> ab_normal_in_abc = a_to_b.cross(abc_normal);
  assert(a_to_c.dot(ab_normal_in_abc) <= 0);
  if (o_to_a.dot(ab_normal_in_abc) <= 0) {
    // ON AB
    simplex.vertices[0] = B;
    simplex.vertices[1] = A;
    simplex.rank = 2;

    // direction is the normal of AB on OAB plane
    const Vector3<T> ab_cross_ao = o_to_a.cross(a_to_b);
    direction = ab_cross_ao.cross(a_to_b);
    direction.normalize();
    assert(o_to_a.dot(direction) <= 0);
    return ProjectionStatus::Continue;
  }

  // On triangle ABC, do not touch the simplex
  const T area = abc_normal.norm();
  if (area < tolerance_ * tolerance_) {
    return ProjectionStatus::Failed;
  }
  const Vector3<T> abc_normal_unit = abc_normal / area;
  const T abc_normal_unit_dot_oa = abc_normal_unit.dot(a);
  if (std::abs(abc_normal_unit_dot_oa) < tolerance_) {
    return ProjectionStatus::Intersect;
  }

  // Use the normal as search direction
  if (abc_normal_unit_dot_oa <= 0) {
    direction = abc_normal_unit;
    return ProjectionStatus::Continue;
  } else {
    direction = -abc_normal_unit;
    return ProjectionStatus::Continue;
  }
}

template <typename T>
typename GJK<T>::ProjectionStatus GJK<T>::simplexProjection4(
    GJKSimplex<T>& simplex, Vector3<T>& direction) const {
  // A is the newly added vertex
  MinkowskiDiffVertex<T> D = simplex.vertices[0];
  MinkowskiDiffVertex<T> C = simplex.vertices[1];
  MinkowskiDiffVertex<T> B = simplex.vertices[2];
  MinkowskiDiffVertex<T> A = simplex.vertices[3];

  const Vector3<T>& a = A.vertex;
  const Vector3<T>& b = B.vertex;
  const Vector3<T>& c = C.vertex;
  const Vector3<T>& d = D.vertex;

  const Vector3<T> ab = b - a;
  const Vector3<T> ac = c - a;
  const Vector3<T> ad = d - a;
  Vector3<T> abc_normal = ab.cross(ac);
  Vector3<T> acd_normal = ac.cross(ad);
  Vector3<T> abd_normal = ab.cross(ad);
  T abc_normal_dot_ad = abc_normal.dot(ad);
  T acd_normal_dot_ab = acd_normal.dot(ab);
  T abd_normal_dot_ac = abd_normal.dot(ac);

  // Check zero
  if (std::abs(abc_normal_dot_ad) <= 0.0) {
    // No volume
    return ProjectionStatus::FailedTetrahedronZeroVolume;
  }

  // Re-orient the normal
  if (abc_normal_dot_ad > 0) {
    abc_normal *= -1;
    abc_normal_dot_ad *= -1;
  }
  if (acd_normal_dot_ab > 0) {
    acd_normal *= -1;
    acd_normal_dot_ab *= -1;
  }
  if (abd_normal_dot_ac > 0) {
    abd_normal *= -1;
    abd_normal_dot_ac *= -1;
  }

  // Check separation
  const Vector3<T>& o_to_a = a;
  const bool d_o_on_same_side_of_abc = o_to_a.dot(abc_normal) > 0;
  const bool c_o_on_same_side_of_abd = o_to_a.dot(abd_normal) > 0;
  const bool b_o_on_same_side_of_acd = o_to_a.dot(acd_normal) > 0;
  if (d_o_on_same_side_of_abc && c_o_on_same_side_of_abd &&
      b_o_on_same_side_of_acd) {
    return ProjectionStatus::Intersect;
  }

  // The libccd implementation: not correct, but let's try
  // Update: surprisingly it works well
  if (!b_o_on_same_side_of_acd) {
    // Remove b
    simplex.vertices[2] = A;
    simplex.rank = 3;
    return simplexProjection3(simplex, direction);
  } else if (!c_o_on_same_side_of_abd) {
    // Remove c
    simplex.vertices[1] = B;
    simplex.vertices[2] = A;
    simplex.rank = 3;
    return simplexProjection3(simplex, direction);
  } else {
    // Remove d
    simplex.vertices[0] = C;
    simplex.vertices[1] = B;
    simplex.vertices[2] = A;
    simplex.rank = 3;
    return simplexProjection3(simplex, direction);
  }
}

}  // namespace cvx_collide
}  // namespace fcl