//
// Created by wei on 23-1-12.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
bool GJK<T>::findMinimumDistancePointsWithSeparatedVertexInit(
    const MinkowskiDiff<T>& shape, GJKSimplex<T>& simplex,
    std::pair<Vector3<T>, Vector3<T>>& p0p1_in_frame0) const {
  // Contains the witness point
  assert(simplex.rank == 1);
  if (simplex.rank != 1) {
    return false;
  }

  // Iteration variable except simplex
  Vector3<T> current_min_distance_point;
  Vector3<T> next_direction;
  MinkowskiDiffVertex<T> new_vertex;
  T bookkeeping_min_distance;

  // Init the current_min_distance_point
  auto init_status =
      computeMinDistanceAndUpdateSimplex(simplex, current_min_distance_point);
  if (init_status != MinDistanceUpdateStatus::OK) return false;

  // Init the distance and direction
  bookkeeping_min_distance = current_min_distance_point.norm();
  if (bookkeeping_min_distance <= tolerance_) {
    return extractSeparationPointNoSubSimplex(shape, simplex, p0p1_in_frame0);
  }
  // The direction is the negative of min_distance_point
  // As it is not norm-zero, we can normalize it
  next_direction = -current_min_distance_point / bookkeeping_min_distance;

  // Start the loop
  std::size_t iteration = 0;
  const T tolerance_squared = tolerance_ * tolerance_;
  while (iteration < max_iterations_) {
    // Update the index
    iteration += 1;

    // Compute the new vertex
    new_vertex = shape.supportVertex(next_direction);

    // Terminate checking on new vertex
    {
      // 1: Should not be very close to the original simplex
      const T delta_on_direction =
          next_direction.dot(new_vertex.vertex - current_min_distance_point);
      if (delta_on_direction < tolerance_) {
        return extractSeparationPointNoSubSimplex(shape, simplex,
                                                  p0p1_in_frame0);
      }

      // 2: Should not be the same as old ones
      for (auto j = 0; j < simplex.n_vertices(); j++) {
        const Vector3<T>& vertex_j = simplex.vertices[j].vertex;
        if ((vertex_j - new_vertex.vertex).squaredNorm() < tolerance_squared) {
          return extractSeparationPointNoSubSimplex(shape, simplex,
                                                    p0p1_in_frame0);
        }
      }
    }

    // Update the simplex
    assert(simplex.rank >= 1 && simplex.rank <= 3);
    simplex.AddVertex(new_vertex);

    // Compute the min-distance point and update the simplex to match that
    auto update_status =
        computeMinDistanceAndUpdateSimplex(simplex, current_min_distance_point);
    assert(simplex.rank >= 1 && simplex.rank <= 3);

    // Check the termination
    if (update_status == MinDistanceUpdateStatus::NoImprovement) {
      return extractSeparationPointNoSubSimplex(shape, simplex, p0p1_in_frame0);
    } else if (update_status == MinDistanceUpdateStatus::OK) {
      // Check if the update converges by numerical
      const T new_min_distance = current_min_distance_point.norm();
      const T distance_improvement =
          bookkeeping_min_distance - new_min_distance;
      if (distance_improvement < tolerance_) {
        return extractSeparationPointNoSubSimplex(shape, simplex,
                                                  p0p1_in_frame0);
      } else if (new_min_distance < tolerance_) {
        return extractSeparationPointNoSubSimplex(shape, simplex,
                                                  p0p1_in_frame0);
      }

      // Continue on the next iteration
      bookkeeping_min_distance = new_min_distance;
      next_direction = -current_min_distance_point / bookkeeping_min_distance;
    } else {
      return false;
    }
  }

  // Reach the iteration limit, claim invalid
  return false;
}

template <typename T>
typename GJK<T>::MinDistanceUpdateStatus
GJK<T>::computeMinDistanceAndUpdateSimplex(
    GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const {
  if (simplex.rank == 1) {
    min_distance_output = simplex.vertices[0].vertex;
    return MinDistanceUpdateStatus::OK;
  } else if (simplex.rank == 2) {
    computeMinDistanceAndUpdateSimplex2(simplex, min_distance_output);
    return MinDistanceUpdateStatus::OK;
  } else if (simplex.rank == 3) {
    computeMinDistanceAndUpdateSimplex3(simplex, min_distance_output);
    return MinDistanceUpdateStatus::OK;
  } else if (simplex.rank == 4) {
    return computeMinDistanceAndUpdateSimplex4(simplex, min_distance_output);
  } else {
    return MinDistanceUpdateStatus::Failed;
  }
}

template <typename T>
void GJK<T>::computeMinDistanceAndUpdateSimplex2(
    GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const {
  assert(simplex.rank == 2);
  // The last one is the new vertex
  const Vector3<T>& s1_new = simplex.vertices[1].vertex;
  const Vector3<T>& s2 = simplex.vertices[0].vertex;

  // Project to s1
  const Vector3<T> s1_to_s2 = s2 - s1_new;
  const T squared_length = s1_to_s2.squaredNorm();
  const T t = -s1_new.dot(s1_to_s2);
  if (t <= 0 || squared_length <= tolerance_ * tolerance_) {
    // Update to point simplex of s1
    min_distance_output = s1_new;
    simplex.vertices[0] = simplex.vertices[1];
    simplex.rank = 1;
  } else if (t >= squared_length) {
    // Update to point simplex of s2
    min_distance_output = s2;
    simplex.rank = 1;
  } else {
    // To segment, the simplex is the same as before
    const T s2_weight = (t / squared_length);
    min_distance_output = s2_weight * s2 + (T(1.0) - s2_weight) * s1_new;
  }
}

template <typename T>
void GJK<T>::computeMinDistanceAndUpdateSimplex3(
    GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const {
  assert(simplex.rank == 3);
  // The last one is the new vertex
  const Vector3<T>& s1_new = simplex.vertices[2].vertex;
  const Vector3<T>& s2 = simplex.vertices[1].vertex;
  const Vector3<T>& s3 = simplex.vertices[0].vertex;

  // The case of s1
  const Vector3<T> s1_to_s2 = s2 - s1_new;
  const Vector3<T> s1_to_s3 = s3 - s1_new;
  const bool s1_separate_s2_o = s1_new.dot(s1_to_s2) >= 0;
  const bool s1_separate_s3_o = s1_new.dot(s1_to_s3) >= 0;
  if (s1_separate_s2_o && s1_separate_s3_o) {
    // Update to point simplex of s1
    min_distance_output = s1_new;
    simplex.vertices[0] = simplex.vertices[2];
    simplex.rank = 1;
    return;
  }

  // The line separation condition needs normal
  const Vector3<T> s1s2s3_plane_normal = s1_to_s2.cross(s1_to_s3);
  const T area_squared = s1s2s3_plane_normal.squaredNorm();
  const bool is_zero_area = area_squared <= T(0.0);

  // We can determine it is segment s1s2
  // s1 here is actually o_to_s1 = s1
  const Vector3<T> s1s2_normal_in_s123 = s1s2s3_plane_normal.cross(s1_to_s2);
  // Old impl:
  // const bool s1s2_separate_s3_o = is_sign_matched(
  //     s1_new.dot(s1s2_normal_in_s123), s1_to_s3.dot(s1s2_normal_in_s123));
  const bool s1s2_separate_s3_o = s1_new.dot(s1s2_normal_in_s123) > 0;
  if ((!s1_separate_s2_o) && s1s2_separate_s3_o) {
    // Keep vertex[2, 1]
    simplex.vertices[0] = simplex.vertices[1];
    simplex.vertices[1] = simplex.vertices[2];
    simplex.rank = 2;
    computeMinDistanceAndUpdateSimplex2(simplex, min_distance_output);
    return;
  }

  // We can determine it is segment s1s3
  const Vector3<T> s1s3_normal_in_s123 = s1s2s3_plane_normal.cross(s1_to_s3);
  // Old impl
  // const bool s1s3_separate_s2_o = is_sign_matched(
  //     s1_new.dot(s1s3_normal_in_s123), s1_to_s2.dot(s1s3_normal_in_s123));
  const bool s1s3_separate_s2_o = s1_new.dot(s1s3_normal_in_s123) < 0;
  if ((!s1_separate_s3_o) && s1s3_separate_s2_o) {
    // Keep vertex[2, 0]
    simplex.vertices[1] = simplex.vertices[2];
    simplex.rank = 2;
    computeMinDistanceAndUpdateSimplex2(simplex, min_distance_output);
    return;
  }

  // Here we need to iterate over all edge
  T min_distance_square = -1;
  Vector3<T> min_distance_point_in_edges, min_distance_point_cache;
  GJKSimplex<T> min_distance_simplex, simplex_cache;
  if (is_zero_area || s1s2_separate_s3_o) {
    // Keep vertex[2, 1]
    simplex_cache = simplex;
    simplex_cache.vertices[0] = simplex_cache.vertices[1];
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.rank = 2;
    computeMinDistanceAndUpdateSimplex2(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_s12 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 || min_distance_s12 < min_distance_square) {
      min_distance_square = min_distance_s12;
      min_distance_simplex = simplex_cache;
      min_distance_point_in_edges = min_distance_point_cache;
    }
  }

  if (is_zero_area || s1s3_separate_s2_o) {
    // Keep vertex[2, 0]
    simplex_cache = simplex;
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.rank = 2;
    computeMinDistanceAndUpdateSimplex2(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_s13 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 || min_distance_s13 < min_distance_square) {
      min_distance_square = min_distance_s13;
      min_distance_simplex = simplex_cache;
      min_distance_point_in_edges = min_distance_point_cache;
    }
  }

  // Last case, s23
  const Vector3<T> s2s3_normal_in_s123 = s1s2s3_plane_normal.cross(s3 - s2);
  // Old impl
  // const bool s2s3_separate_s1_o = is_sign_matched(
  //    - s2.dot(s2s3_normal_in_s123), s1_to_s2.dot(s2s3_normal_in_s123));
  const bool s2s3_separate_s1_o = s2.dot(s2s3_normal_in_s123) > 0;
  if (is_zero_area || s2s3_separate_s1_o) {
    // Keep vertex[1, 0]
    simplex_cache = simplex;
    simplex_cache.rank = 2;
    computeMinDistanceAndUpdateSimplex2(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_s23 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 || min_distance_s23 < min_distance_square) {
      min_distance_square = min_distance_s23;
      min_distance_simplex = simplex_cache;
      min_distance_point_in_edges = min_distance_point_cache;
    }
  }

  // Now, do triangle s123 if necessary
  // If min_distance_square is updated, then there must be one edge s_ij that
  // separates s_k and o. In this situation, o can not be projected to the
  // internal of s_ijk (and we do NOT need to check that).
  // On the contrary, if min_distance_square is not updated, then o must be
  // projected into the internal of s_ijk (s_123).
  if (min_distance_square < 0) {
    assert(!is_zero_area);
    const Vector3<T>& n = s1s2s3_plane_normal;
    const T d = s1_new.dot(n);
    const Vector3<T> o_to_project = n * (d / area_squared);
    min_distance_output = o_to_project;
  } else {
    // Just use the smallest distance one
    // Important: we cannot return NoImprovement as this method is used as a
    //            subroutine for simplex4.
    //            If it is not a subroutine, then edge s23 achieves min-distance
    //            implies no improvement.
    simplex = min_distance_simplex;
    min_distance_output = min_distance_point_in_edges;
  }
}

template <typename T>
typename GJK<T>::MinDistanceUpdateStatus
GJK<T>::computeMinDistanceAndUpdateSimplex4(
    GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const {
  assert(simplex.rank == 4);
  // The last one is the new vertex
  // const Vector3<T>& s1_new = simplex.vertices[3].vertex;
  // const Vector3<T>& s2 = simplex.vertices[2].vertex;
  // const Vector3<T>& s3 = simplex.vertices[1].vertex;
  // const Vector3<T>& s4 = simplex.vertices[0].vertex;

  // As we know it is seperated, the min-distance must be one of its triangle
  // In particular s123, s124, s134
  // The cache output
  T min_distance_square = -1;
  Vector3<T> min_distance_point_in_faces, min_distance_point_cache;
  GJKSimplex<T> min_distance_simplex, simplex_cache;

  // s123
  {
    simplex_cache = simplex;
    // Keep vertex[3, 2, 1]
    simplex_cache.vertices[0] = simplex_cache.vertices[1];
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.vertices[2] = simplex_cache.vertices[3];
    simplex_cache.rank = 3;

    // We need to run weight checking as in sub-simplex, the min-distance
    // triangle may NOT contain s1
    computeMinDistanceAndUpdateSimplex3(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_square_s123 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 ||
        min_distance_square_s123 < min_distance_square) {
      min_distance_square = min_distance_square_s123;
      min_distance_point_in_faces = min_distance_point_cache;
      min_distance_simplex = simplex_cache;
    }
  }

  // s124
  {
    simplex_cache = simplex;
    // Keep vertex[3, 2, 0]
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.vertices[2] = simplex_cache.vertices[3];
    simplex_cache.rank = 3;
    computeMinDistanceAndUpdateSimplex3(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_square_s124 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 ||
        min_distance_square_s124 < min_distance_square) {
      min_distance_square = min_distance_square_s124;
      min_distance_point_in_faces = min_distance_point_cache;
      min_distance_simplex = simplex_cache;
    }
  }

  // s134
  {
    simplex_cache = simplex;
    // Keep vertex[3, 1, 0]
    simplex_cache.vertices[2] = simplex_cache.vertices[3];
    simplex_cache.rank = 3;
    computeMinDistanceAndUpdateSimplex3(simplex_cache,
                                        min_distance_point_cache);
    const T min_distance_square_s134 = min_distance_point_cache.squaredNorm();
    if (min_distance_square < 0 ||
        min_distance_square_s134 < min_distance_square) {
      min_distance_square = min_distance_square_s134;
      min_distance_point_in_faces = min_distance_point_cache;
      min_distance_simplex = simplex_cache;
    }
  }

  // Improvement
  if (min_distance_square < 0) {
    return MinDistanceUpdateStatus::NoImprovement;
  } else {
    simplex = min_distance_simplex;
    min_distance_output = min_distance_point_in_faces;
    return MinDistanceUpdateStatus::OK;
  }
}

template <typename T>
bool GJK<T>::extractSeparationPointNoSubSimplex(
    const MinkowskiDiff<T>& shape, const GJKSimplex<T>& simplex,
    std::pair<Vector3<T>, Vector3<T>>& output) const {
  // The ordering (old/new) of vertices does NOT matter in this method
  if ((simplex.rank == 4) || (!simplex.is_valid())) {
    // Must be valid simplex, and the size cannot be 4
    return false;
  }

  // Simple case
  if (simplex.rank == 1) {
    const auto& v = simplex.vertices[0];
    output.first = shape.support0(v.direction);
    output.second = shape.support1(-v.direction);
    return true;
  } else if (simplex.rank == 2) {
    const auto& v1 = simplex.vertices[0];
    const auto& v2 = simplex.vertices[1];
    const Vector3<T>& s1_new = v1.vertex;
    const Vector3<T>& s2 = v2.vertex;
    const Vector3<T> s1_to_s2 = s2 - s1_new;
    const T squared_length = s1_to_s2.squaredNorm();

    // This should not happen if we reach here, but let's handle it
    if (squared_length <= 0.0) {
      const auto& v = simplex.vertices[0];
      output.first = shape.support0(v.direction);
      output.second = shape.support1(-v.direction);
      return true;
    }

    // Now we can do division
    const T t = -s1_new.dot(s1_to_s2);
    const T s2_weight = t / squared_length;

    // Should be in [0, 1]
    if (s2_weight > 1 + barycentric_weight_tolerance) {
      // Should be s2, write output but return false
      output.first = shape.support0(v2.direction);
      output.second = shape.support1(-v2.direction);
      return false;
    } else if (s2_weight < -barycentric_weight_tolerance) {
      // Should be s1
      output.first = shape.support0(v1.direction);
      output.second = shape.support1(-v1.direction);
      return false;
    } else {
      // Assign the output
      const T s1_weight = T(1.0) - s2_weight;
      output.first = shape.support0(v1.direction) * s1_weight +
                     shape.support0(v2.direction) * s2_weight;
      output.second = shape.support1(-v1.direction) * s1_weight +
                      shape.support1(-v2.direction) * s2_weight;
      return true;
    }
  } else {
    assert(simplex.rank == 3);
    const auto& v1 = simplex.vertices[0];
    const auto& v2 = simplex.vertices[1];
    const auto& v3 = simplex.vertices[2];
    const Vector3<T>& s1 = v1.vertex;
    const Vector3<T>& s2 = v2.vertex;
    const Vector3<T>& s3 = v3.vertex;

    // Compute the weight
    const Vector3<T> s1_to_s2 = s2 - s1;
    const Vector3<T> s1_to_s3 = s3 - s1;
    const Vector3<T> s1s2s3_plane_normal = s1_to_s2.cross(s1_to_s3);
    const T area_squared = s1s2s3_plane_normal.squaredNorm();
    if (area_squared <= 0.0) return false;  // Should not happen if we are here
    const Vector3<T>& n = s1s2s3_plane_normal;
    const T d = s1.dot(n);
    const Vector3<T> o_to_project = n * (d / area_squared);
    const T area = std::sqrt(area_squared);
    const T s2_weight = (s1_to_s3.cross(s1 - o_to_project)).norm() / area;
    const T s3_weight = (s1_to_s2.cross(s1 - o_to_project)).norm() / area;
    const T s1_weight = T(1.0) - s2_weight - s3_weight;

    // All should be in [0, 1]
    if ((s1_weight < -barycentric_weight_tolerance) ||
        (s2_weight < -barycentric_weight_tolerance) ||
        (s3_weight < -barycentric_weight_tolerance)) {
      return false;
    }

    // Assign the output
    output.first = shape.support0(v1.direction) * s1_weight +
                   shape.support0(v2.direction) * s2_weight +
                   shape.support0(v3.direction) * s3_weight;
    output.second = shape.support1(-v1.direction) * s1_weight +
                    shape.support1(-v2.direction) * s2_weight +
                    shape.support1(-v3.direction) * s3_weight;
    return true;
  }
}

template <typename T>
bool GJK<T>::extractSeparationPointTrySubSimplex(
    const MinkowskiDiff<T>& shape, const GJKSimplex<T>& simplex,
    std::pair<Vector3<T>, Vector3<T>>& output) const {
  if (simplex.rank >= 4) {
    return false;
  } else if (simplex.rank <= 2) {
    extractSeparationPointNoSubSimplex(shape, simplex, output);
    // Always OK here, despite the subtle when rank==2
    return true;
  }

  // Must be 3-simplex
  assert(simplex.rank == 3);
  T min_distance_squared = -1;
  GJKSimplex<T> simplex_cache;
  std::pair<Vector3<T>, Vector3<T>> min_distance_output, output_cache;

  // Triangle 123
  {
    simplex_cache = simplex;
    auto ok =
        extractSeparationPointNoSubSimplex(shape, simplex_cache, output_cache);
    if (ok) {
      // Must be on the triangle
      const auto this_distance_square =
          (output_cache.first - output_cache.second).squaredNorm();
      if (min_distance_squared < 0 ||
          (this_distance_square < min_distance_squared)) {
        min_distance_squared = this_distance_square;
        min_distance_output = output_cache;
      }
    }
  }

  // Segment 12
  {
    simplex_cache = simplex;
    simplex_cache.rank = 2;
    extractSeparationPointNoSubSimplex(shape, simplex_cache, output_cache);
    // Must be on valid
    const auto this_distance_square =
        (output_cache.first - output_cache.second).squaredNorm();
    if (min_distance_squared < 0 ||
        (this_distance_square < min_distance_squared)) {
      min_distance_squared = this_distance_square;
      min_distance_output = output_cache;
    }
  }

  // Segment 13
  {
    simplex_cache = simplex;
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.rank = 2;
    extractSeparationPointNoSubSimplex(shape, simplex_cache, output_cache);
    // Must be on valid
    const auto this_distance_square =
        (output_cache.first - output_cache.second).squaredNorm();
    if (min_distance_squared < 0 ||
        (this_distance_square < min_distance_squared)) {
      min_distance_squared = this_distance_square;
      min_distance_output = output_cache;
    }
  }

  // Segment 23
  {
    simplex_cache = simplex;
    simplex_cache.vertices[0] = simplex_cache.vertices[1];
    simplex_cache.vertices[1] = simplex_cache.vertices[2];
    simplex_cache.rank = 2;
    extractSeparationPointNoSubSimplex(shape, simplex_cache, output_cache);
    // Must be on valid
    const auto this_distance_square =
        (output_cache.first - output_cache.second).squaredNorm();
    if (min_distance_squared < 0 ||
        (this_distance_square < min_distance_squared)) {
      min_distance_squared = this_distance_square;
      min_distance_output = output_cache;
    }
  }

  // Assign the output
  assert(min_distance_squared >= 0);
  output = min_distance_output;
  return true;
}

}  // namespace cvx_collide
}  // namespace fcl