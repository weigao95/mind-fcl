//
// Created by wei on 22-6-5.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
typename EPA<T>::FindNextDirectionStatus EPA<T>::findNextSupportDirection(
    const Polytope<T>& polytope, const MinkowskiDiff<T>& shape,
    PolytopeElementBase* nearest_feature, T tolerance, Vector3<T>& next_d,
    Vector3<T>& next_v, PolytopeFace<T>** start_face) const {
  auto feature_type = nearest_feature->type;
  // Usually we select the min-distance point on the feature
  // But when the feature is too close, we might choose the outer normal
  // This happens when GJK terminate upon simplex closet to origin
  constexpr T use_outer_normal_if_nearest_smaller_than = T(1e-3);
  constexpr T outer_normal_threshold =
      use_outer_normal_if_nearest_smaller_than *
      use_outer_normal_if_nearest_smaller_than;

  // Check face
  auto process_face_candidate =
      [&](PolytopeFace<T>* candidate_face,
          bool try_min_distance_witness) -> FindNextDirectionStatus {
    // Compute face normal
    Vector3<T> face_normal;
    auto normal_ok =
        polytope.ComputeFaceNormalPointingOutward(candidate_face, face_normal);
    if (!normal_ok) return FindNextDirectionStatus::Failed;

    // Get a point on the face
    const Vector3<T>& point_on_face =
        candidate_face->edges_of_face[0]->vertices_of_edge[0]->vertex_location;

    // First try min-distance
    const auto& min_distance = candidate_face->min_distance;
    if (try_min_distance_witness &&
        (min_distance.min_distance_square > outer_normal_threshold)) {
      // Case 1: directly search on the witness point
      next_d = min_distance.witness;
      next_d.normalize();

      // Compute the point
      next_v = shape.support(next_d);
      const T next_v_delta = face_normal.dot(next_v - point_on_face);

      // Check the outside is OK
      if (next_v_delta >= T(0.0)) {
        // OK
        *start_face = candidate_face;
        return FindNextDirectionStatus::OK;
      }
    }

    // We should try the normal
    next_d = face_normal;
    next_v = shape.support(next_d);

    // This is the distance on the normal direction
    const T next_v_delta = face_normal.dot(next_v - point_on_face);

    // Check the outside is OK
    if (next_v_delta > tolerance) {
      // OK
      *start_face = candidate_face;
      return FindNextDirectionStatus::OK;
    } else {
      // We can determine converge here
      return FindNextDirectionStatus::Converge;
    }
  };

  // According to the nearest feature type
  if (feature_type == PolytopeElementType::Face) {
    auto* face = static_cast<PolytopeFace<T>*>(nearest_feature);
    return process_face_candidate(face, true);
  } else if (feature_type == PolytopeElementType::Edge) {
    auto* edge = static_cast<PolytopeEdge<T>*>(nearest_feature);
    auto* face_0 = edge->faces_of_edge[0];
    auto* face_1 = edge->faces_of_edge[1];
    const auto& min_distance = edge->min_distance;
    if (min_distance.min_distance_square > outer_normal_threshold) {
      // Case 1: directly search on the witness point
      next_d = min_distance.witness;
      next_d.normalize();

      // Compute the point
      next_v = shape.support(next_d);

      // Check the outside is OK for face_0
      if (polytope.IsPointOutsidePolytopeFace(face_0, next_v)) {
        *start_face = face_0;
        return FindNextDirectionStatus::OK;
      }

      // Now do face_1
      if (polytope.IsPointOutsidePolytopeFace(face_1, next_v)) {
        *start_face = face_1;
        return FindNextDirectionStatus::OK;
      }
    }

    // For faces on the edge
    auto face_0_status = process_face_candidate(face_0, false);
    if (face_0_status == FindNextDirectionStatus::OK) return face_0_status;
    return process_face_candidate(face_1, false);
  } else {
    return FindNextDirectionStatus::Failed;
  }
}

template <typename T>
EPA_Status EPA<T>::evaluateFromUnInitializedPolytope(
    const GJKSimplex<T>& simplex, Polytope<T>& polytope,
    const MinkowskiDiff<T>& shape, T* depth_if_penetration,
    Vector3<T>* p_on_shape_0, Vector3<T>* p_on_shape_1) const {
  auto simple2polytope_status = simplexToPolytope(
      simplex, shape, polytope, p_on_shape_0, p_on_shape_1, tolerance_);

  // Depends on the status
  if (simple2polytope_status == Simplex2PolytopeStatus::Failed) {
    return EPA_Status::Failed;
  } else if (simple2polytope_status == Simplex2PolytopeStatus::Touching) {
    if (depth_if_penetration) *depth_if_penetration = 0.0;
    return EPA_Status::Touching;
  }

  // Forward to impl
  return evaluateFromInitializedPolytope(polytope, shape, depth_if_penetration,
                                         p_on_shape_0, p_on_shape_1);
}

template <typename T>
EPA_Status EPA<T>::evaluateFromInitializedPolytope(
    Polytope<T>& polytope, const MinkowskiDiff<T>& shape,
    T* depth_if_penetration, Vector3<T>* p_on_shape_0,
    Vector3<T>* p_on_shape_1) const {
  // The iterations
  std::size_t iteration = 0;
  RawFeature min_distance_raw;
  while (true) {
    // Find the nearest feature
    auto* min_distance_feature = polytope.ComputeMinDistanceToOrigin();

    // Should be valid
    if (min_distance_feature == nullptr) return EPA_Status::Failed;

    // Should be edge/face
    auto feature_type = min_distance_feature->type;
    if (feature_type == PolytopeElementType::Vertex) {
      // Usually it should not be a vertex, but there might exist very
      // coincidental cases?
      auto* vertex = static_cast<PolytopeVertex<T>*>(min_distance_feature);
      const auto& min_distance = vertex->min_distance;
      constexpr T ratio = T(1) - T(1e-3);
      auto* other_feature = polytope.ComputeMinDistanceToOrigin(true);
      if (other_feature == nullptr) return EPA_Status::Failed;
      if (other_feature->type == PolytopeElementType::Edge) {
        auto* edge = static_cast<PolytopeEdge<T>*>(other_feature);
        if (ratio * (edge->min_distance.min_distance_square) <
            min_distance.min_distance_square) {
          min_distance_feature = edge;
          feature_type = min_distance_feature->type;
        }
      } else if (other_feature->type == PolytopeElementType::Face) {
        auto* face = static_cast<PolytopeFace<T>*>(other_feature);
        if (ratio * (face->min_distance.min_distance_square) <
            min_distance.min_distance_square) {
          min_distance_feature = face;
          feature_type = min_distance_feature->type;
        }
      }
    }

    // Should be edge/face here
    if ((feature_type != PolytopeElementType::Edge) &&
        (feature_type != PolytopeElementType::Face)) {
      return EPA_Status::Failed;
    }

    // To raw feature
    auto ok = transformToRawFeature(min_distance_feature, min_distance_raw);
    assert(ok);
    (void)(ok);  // disable warning

    // The direction
    Vector3<T> next_direction, next_v;
    PolytopeFace<T>* start_face;
    auto direction_status = findNextSupportDirection(
        polytope, shape, min_distance_feature, tolerance_, next_direction,
        next_v, &start_face);

    // Check the status of next direction
    if (direction_status == FindNextDirectionStatus::Failed) {
      return EPA_Status::Failed;
    } else if (direction_status == FindNextDirectionStatus::Converge) {
      assignPenetrationPair(shape, {next_v, next_direction}, min_distance_raw,
                            depth_if_penetration, p_on_shape_0, p_on_shape_1);
      return EPA_Status::OK;
    }

    // Here we should get FindNextDirectionStatus::OK
    if (start_face == nullptr) return EPA_Status::Failed;
    assert(start_face != nullptr);
    MinkowskiDiffVertex<T> next_vertex{next_v, next_direction};

    // Check termination
    auto converge = checkTerminateCondition(shape, min_distance_raw,
                                            next_direction, next_v, tolerance_);
    if (converge) {
      assignPenetrationPair(shape, next_vertex, min_distance_raw,
                            depth_if_penetration, p_on_shape_0, p_on_shape_1);
      return EPA_Status::OK;
    }

    // Expand the polytope
    auto expand_status = polytope.ExpandPolytope(next_vertex, start_face);
    if (expand_status == Polytope<T>::PolytopeExpansionStatus::Failed) {
      return EPA_Status::Failed;
    } else if (expand_status ==
               Polytope<T>::PolytopeExpansionStatus::MallocFailed) {
      assignPenetrationPair(shape, next_vertex, min_distance_raw,
                            depth_if_penetration, p_on_shape_0, p_on_shape_1);
      return EPA_Status::MallocFailed;
    }

    iteration += 1;
    if (iteration >= max_iterations_) {
      assignPenetrationPair(shape, next_vertex, min_distance_raw,
                            depth_if_penetration, p_on_shape_0, p_on_shape_1);
      return EPA_Status::IterationLimit;
    }
  }
}

template <typename T>
EPA_Status EPA<T>::Evaluate(const GJKSimplex<T>& simplex,
                            const MinkowskiDiff<T>& shape,
                            T* depth_if_penetration, Vector3<T>* p_on_shape_0,
                            Vector3<T>* p_on_shape_1,
                            Polytope<T>* external_polytope) const {
  // Make the polytopes
  if (external_polytope == nullptr) {
    Polytope<T> polytope(max_n_faces_);
    return evaluateFromUnInitializedPolytope(simplex, polytope, shape,
                                             depth_if_penetration, p_on_shape_0,
                                             p_on_shape_1);
  } else {
    external_polytope->Reset(max_n_faces_);
    return evaluateFromUnInitializedPolytope(simplex, *external_polytope, shape,
                                             depth_if_penetration, p_on_shape_0,
                                             p_on_shape_1);
  }
}

template <typename T>
EPA_Status EPA<T>::Test_EvaluateFromInitializedPolytope(
    Polytope<T>& polytope, const MinkowskiDiff<T>& shape,
    T* depth_if_penetration, Vector3<T>* p_on_shape_0,
    Vector3<T>* p_on_shape_1) const {
  return evaluateFromInitializedPolytope(polytope, shape, depth_if_penetration,
                                         p_on_shape_0, p_on_shape_1);
}

template <typename T>
bool EPA<T>::checkTerminateCondition(const MinkowskiDiff<T>& shape,
                                     const RawFeature& nearest_feature,
                                     const Vector3<T>& d,
                                     const Vector3<T>& new_v,
                                     T tolerance) const {
  (void)shape;
  (void)d;
  auto feature_type = nearest_feature.type;
  assert(feature_type == PolytopeElementType::Edge ||
         feature_type == PolytopeElementType::Face);

  // Check the changing distance
  T delta_distance_squared = T(0.0);
  if (feature_type == PolytopeElementType::Edge) {
    const Vector3<T>& a = nearest_feature.a_vertex.vertex;
    const Vector3<T>& b = nearest_feature.b_vertex.vertex;
    bool in_segment{false};
    pointToSegmentSquaredDistance<T>(new_v, a, b, in_segment,
                                     delta_distance_squared, nullptr);
  } else {
    const Vector3<T> a = nearest_feature.a_vertex.vertex;
    const Vector3<T> b = nearest_feature.b_vertex.vertex;
    const Vector3<T> c = nearest_feature.c_vertex.vertex;
    const T point2plane = pointToPlaneDistance<T>(new_v, a, b, c);
    delta_distance_squared = point2plane * point2plane;
  }

  // If ok
  bool converge = (delta_distance_squared < tolerance * tolerance);
  return converge;
}

template <typename T>
bool EPA<T>::transformToRawFeature(const PolytopeElementBase* feature,
                                   RawFeature& raw_feature) {
  if (feature->type == PolytopeElementType::Vertex) {
    const auto* vertex = static_cast<const PolytopeVertex<T>*>(feature);
    raw_feature.type = PolytopeElementType::Vertex;
    raw_feature.min_distance = vertex->min_distance;

    // For vertex
    raw_feature.a_vertex.vertex = vertex->vertex_location;
    raw_feature.a_vertex.direction = vertex->minkowski_diff_direction;
    return true;
  } else if (feature->type == PolytopeElementType::Edge) {
    auto* edge = static_cast<const PolytopeEdge<T>*>(feature);
    raw_feature.type = PolytopeElementType::Edge;
    raw_feature.min_distance = edge->min_distance;

    // For edge
    const Vector3<T>& a = edge->vertices_of_edge[0]->vertex_location;
    const Vector3<T>& b = edge->vertices_of_edge[1]->vertex_location;
    const Vector3<T>& d_a = edge->vertices_of_edge[0]->minkowski_diff_direction;
    const Vector3<T>& d_b = edge->vertices_of_edge[1]->minkowski_diff_direction;
    raw_feature.a_vertex.vertex = a;
    raw_feature.a_vertex.direction = d_a;
    raw_feature.b_vertex.vertex = b;
    raw_feature.b_vertex.direction = d_b;
    return true;
  } else if (feature->type == PolytopeElementType::Face) {
    auto* face = static_cast<const PolytopeFace<T>*>(feature);
    raw_feature.type = PolytopeElementType::Face;
    raw_feature.min_distance = face->min_distance;
    return face->GetVertices(raw_feature.a_vertex, raw_feature.b_vertex,
                             raw_feature.c_vertex);
  }

  // Not suitable
  return false;
}

template <typename T>
void EPA<T>::assignPenetrationPair(
    const MinkowskiDiff<T>& shape,
    const MinkowskiDiffVertex<T>& candidate_next_v,
    const RawFeature& nearest_feature, T* depth_if_penetration,
    Vector3<T>* p_on_shape_0, Vector3<T>* p_on_shape_1) const {
  // Depends on the type of feature
  PolytopeElementType feature_type = nearest_feature.type;
  if (feature_type == PolytopeElementType::Vertex) {
    const auto& d = nearest_feature.a_vertex.direction;
    if (p_on_shape_0) *p_on_shape_0 = shape.support0(d);
    if (p_on_shape_1) *p_on_shape_1 = shape.support1(-d);
    if (depth_if_penetration)
      *depth_if_penetration = nearest_feature.a_vertex.vertex.norm();
    return;
  } else if (feature_type == PolytopeElementType::Edge) {
    auto ok = assignPenetrationPairFromSegment(shape, nearest_feature,
                                               depth_if_penetration,
                                               p_on_shape_0, p_on_shape_1);
    if (ok) return;
    // else, fall-back
  } else if (feature_type == PolytopeElementType::Face) {
    auto ok = assignPenetrationPairFromFace(shape, nearest_feature,
                                            depth_if_penetration, p_on_shape_0,
                                            p_on_shape_1);
    if (ok) return;
  }

  //  The fall-back implementation
  Vector3<T> p0 = shape.support0(candidate_next_v.direction);
  Vector3<T> p1 = shape.support1(-candidate_next_v.direction);
  if (p_on_shape_0) *p_on_shape_0 = p0;
  if (p_on_shape_1) *p_on_shape_1 = p1;
  if (depth_if_penetration) *depth_if_penetration = (p0 - p1).norm();
}

template <typename T>
bool EPA<T>::assignPenetrationPairFromSegment(
    const MinkowskiDiff<T>& shape, const RawFeature& nearest_edge_feature,
    T* depth_if_penetration, Vector3<T>* p_on_shape_0,
    Vector3<T>* p_on_shape_1) const {
  // The distance is rather simple
  const auto& min_distance_record = nearest_edge_feature.min_distance;
  if (depth_if_penetration) {
    *depth_if_penetration = std::sqrt(min_distance_record.min_distance_square);
  }

  // Use interpolation to compute the point on shape 0/1
  // p = a + s * a_to_b, where 0 <= s <= 1
  const Vector3<T>& p = min_distance_record.witness;
  const Vector3<T>& a = nearest_edge_feature.a_vertex.vertex;
  const Vector3<T>& b = nearest_edge_feature.b_vertex.vertex;
  const Vector3<T>& d_a = nearest_edge_feature.a_vertex.direction;
  const Vector3<T>& d_b = nearest_edge_feature.b_vertex.direction;

  // Compute the index over xyz that a_to_b[i] has largest length
  const Vector3<T> a_to_b = b - a;
  int max_over_xyz_of_ab = 0;
  T max_coordinate_diff_ab = std::abs(a_to_b[0]);
  if (std::abs(a_to_b[1]) > max_coordinate_diff_ab) {
    max_over_xyz_of_ab = 1;
    max_coordinate_diff_ab = std::abs(a_to_b[1]);
  }
  if (std::abs(a_to_b[2]) > max_coordinate_diff_ab) {
    max_over_xyz_of_ab = 2;
    max_coordinate_diff_ab = std::abs(a_to_b[2]);
  }

  // Compute s
  if (max_coordinate_diff_ab <= 1e-10) {
    // Use a?
    if (p_on_shape_0) *p_on_shape_0 = shape.support0(d_a);
    if (p_on_shape_1) *p_on_shape_1 = shape.support1(-d_a);
    return true;
  }
  const T s = (p[max_over_xyz_of_ab] - a[max_over_xyz_of_ab]) /
              (b[max_over_xyz_of_ab] - a[max_over_xyz_of_ab]);
  if (s < 0 || s > 1.0) return false;
  const T a_weight = T(1.0) - s;
  const T b_weight = s;

  // Assign the result
  if (p_on_shape_0) {
    const Vector3<T> v_a_0 = shape.support0(d_a);
    const Vector3<T> v_b_0 = shape.support0(d_b);
    *p_on_shape_0 = a_weight * v_a_0 + b_weight * v_b_0;
  }
  if (p_on_shape_1) {
    const Vector3<T> v_a_1 = shape.support1(-d_a);
    const Vector3<T> v_b_1 = shape.support1(-d_b);
    *p_on_shape_1 = a_weight * v_a_1 + b_weight * v_b_1;
  }
  return true;
}

template <typename T>
bool EPA<T>::assignPenetrationPairFromFace(
    const MinkowskiDiff<T>& shape, const RawFeature& nearest_face_feature,
    T* depth_if_penetration, Vector3<T>* p_on_shape_0,
    Vector3<T>* p_on_shape_1) const {
  // The distance is rather simple
  const auto& min_distance_record = nearest_face_feature.min_distance;
  if (depth_if_penetration) {
    *depth_if_penetration = std::sqrt(min_distance_record.min_distance_square);
  }

  // Get a, b, c
  const MinkowskiDiffVertex<T>& a_vertex = nearest_face_feature.a_vertex;
  const MinkowskiDiffVertex<T>& b_vertex = nearest_face_feature.b_vertex;
  const MinkowskiDiffVertex<T>& c_vertex = nearest_face_feature.c_vertex;

  // Compute the weight
  const Vector3<T>& a = a_vertex.vertex;
  const Vector3<T>& b = b_vertex.vertex;
  const Vector3<T>& c = c_vertex.vertex;
  const Vector3<T>& d_a = a_vertex.direction;
  const Vector3<T>& d_b = b_vertex.direction;
  const Vector3<T>& d_c = c_vertex.direction;

  // Compute the interpolation weight
  std::array<T, 3> parameterization;
  {
    const Vector3<T>& p = min_distance_record.witness;
    const Vector3<T> dl[] = {a - b, b - c, c - a};
    const Vector3<T> n = dl[0].cross(dl[1]);
    const T n_squared_norm = n.squaredNorm();
    const T area = std::sqrt(n_squared_norm);

    // We don't need to check p on abc
    parameterization[0] = dl[1].cross(b - p).norm() / area;
    parameterization[1] = dl[2].cross(c - p).norm() / area;
    auto parameterization_2_for_checking = dl[0].cross(a - p).norm() / area;
    parameterization[2] = T(1.0) - parameterization[0] - parameterization[1];

    // Check the validity on parameterization 2
    const T ratio_tolerance = 0.01;
    if (std::abs(parameterization_2_for_checking - parameterization[2]) >
        ratio_tolerance) {
      return false;
    }
  }

  // Now do assign
  if (p_on_shape_0) {
    const Vector3<T> v_a_0 = shape.support0(d_a);
    const Vector3<T> v_b_0 = shape.support0(d_b);
    const Vector3<T> v_c_0 = shape.support0(d_c);
    *p_on_shape_0 = parameterization[0] * v_a_0 + parameterization[1] * v_b_0 +
                    parameterization[2] * v_c_0;
  }
  if (p_on_shape_1) {
    const Vector3<T> v_a_1 = shape.support1(-d_a);
    const Vector3<T> v_b_1 = shape.support1(-d_b);
    const Vector3<T> v_c_1 = shape.support1(-d_c);
    *p_on_shape_1 = parameterization[0] * v_a_1 + parameterization[1] * v_b_1 +
                    parameterization[2] * v_c_1;
  }
  return true;
}

}  // namespace cvx_collide
}  // namespace fcl