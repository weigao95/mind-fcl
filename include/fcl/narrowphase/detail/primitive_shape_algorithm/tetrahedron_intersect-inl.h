#pragma once

namespace fcl {

namespace detail {
constexpr int kVertexOfTetrahedron = 4;
constexpr int kVertexOfTriangle = 3;

template <typename S>
void calcTetrahedronNormalHeight(
    const std::array<Vector3<S>, kVertexOfTetrahedron>& vertices,
    std::array<Vector3<S>, kVertexOfTetrahedron>& normals,
    std::array<S, kVertexOfTetrahedron>& heights) {
  for (int i = 0; i < kVertexOfTetrahedron; i++) {
    const Vector3<S>& v = vertices[i];
    const Vector3<S>& a = vertices[(i + 1) % kVertexOfTetrahedron];
    const Vector3<S>& b = vertices[(i + 2) % kVertexOfTetrahedron];
    const Vector3<S>& c = vertices[(i + 3) % kVertexOfTetrahedron];

    // Calculate edge vectors and face normal
    const Vector3<S> edge_ab = a - b;
    const Vector3<S> edge_bc = b - c;
    normals[i] = edge_ab.cross(edge_bc);

    // Make sure the normals are facing outwards
    const Vector3<S> edge_av = a - v;
    heights[i] = normals[i].dot(edge_av);

    const bool sign = heights[i] > 0;
    if (sign) {
      heights[i] = -heights[i];
    } else {
      normals[i] = -normals[i];
    }
  }
}

template <typename S, int n_points>
bool checkPointDistanceToFaceAllOutsideBounds(
    const std::array<Vector3<S>, n_points>& points,
    const Vector3<S>& face_normal, const Vector3<S>& point_on_face,
    S lower_bound, S upper_bound) {
  static_assert(n_points >= 1, "points cannot be empty");
  bool at_least_one_smaller_than_ub = false;
  bool at_least_one_larger_than_lb = false;
  for (int i = 0; i < n_points; i++) {
    const auto& edge = points[i] - point_on_face;
    const auto dot = edge.dot(face_normal);
    if (dot < upper_bound) at_least_one_smaller_than_ub = true;
    if (dot > lower_bound) at_least_one_larger_than_lb = true;
    if (at_least_one_larger_than_lb && at_least_one_smaller_than_ub)
      return false;
  }
  assert((!at_least_one_smaller_than_ub) || (!at_least_one_larger_than_lb));
  return true;
}

template <typename S, int n_p1, int n_p2>
bool isSeparatedOnAxis(const std::array<Vector3<S>, n_p1>& p1,
                       const std::array<Vector3<S>, n_p2>& p2,
                       const Vector3<S>& axis) {
  S shape1_lb = std::numeric_limits<S>::infinity();
  S shape1_ub = -std::numeric_limits<S>::infinity();
  for (const auto& v_i : p1) {
    const auto distance_to_O_i = v_i.dot(axis);
    shape1_lb = std::min(shape1_lb, distance_to_O_i);
    shape1_ub = std::max(shape1_ub, distance_to_O_i);
  }

  S shape2_lb = std::numeric_limits<S>::infinity();
  S shape2_ub = -std::numeric_limits<S>::infinity();
  for (const auto& v_i : p2) {
    const auto distance_to_O_i = v_i.dot(axis);
    shape2_lb = std::min(shape2_lb, distance_to_O_i);
    shape2_ub = std::max(shape2_ub, distance_to_O_i);
  }

  const bool separated = shape1_lb > shape2_ub || shape2_lb > shape1_ub;
  return separated;
}

template <typename S>
std::array<Vector3<S>, kVertexOfTetrahedron> getTetrahedronVertices(
    const Tetrahedron<S>& s, const Transform3<S>& tf) {
  std::array<Vector3<S>, kVertexOfTetrahedron> result;
  for (int i = 0; i < kVertexOfTetrahedron; i++) result[i] = tf * s.vertices[i];
  return result;
}

template <typename S>
bool tetrahedronTrahedronIntersect(const Tetrahedron<S>& s1,
                                   const Transform3<S>& tf1,
                                   const Tetrahedron<S>& s2,
                                   const Transform3<S>& tf2) {
  const std::array<Vector3<S>, kVertexOfTetrahedron> vertices1 =
      getTetrahedronVertices(s1, tf1);
  const std::array<Vector3<S>, kVertexOfTetrahedron> vertices2 =
      getTetrahedronVertices(s2, tf2);

  std::array<Vector3<S>, kVertexOfTetrahedron> normals1;
  std::array<S, kVertexOfTetrahedron> height1;
  calcTetrahedronNormalHeight(vertices1, normals1, height1);
  for (int i = 0; i < kVertexOfTetrahedron; i++) {
    const auto& normal_i = normals1[i];
    const auto& point_in_face_i = vertices1[(i + 1) % kVertexOfTetrahedron];

    const bool all_out_bound_i = checkPointDistanceToFaceAllOutsideBounds<S, 4>(
        vertices2, normal_i, point_in_face_i, height1[i], 0);
    if (all_out_bound_i) return false;
  }

  std::array<Vector3<S>, kVertexOfTetrahedron> normals2;
  std::array<S, kVertexOfTetrahedron> height2;
  calcTetrahedronNormalHeight(vertices2, normals2, height2);
  for (int i = 0; i < kVertexOfTetrahedron; i++) {
    const auto& normal_i = normals2[i];
    const auto& point_in_face_i = vertices2[(i + 1) % kVertexOfTetrahedron];

    const bool all_out_bound_i = checkPointDistanceToFaceAllOutsideBounds<S, 4>(
        vertices1, normal_i, point_in_face_i, height2[i], 0);
    if (all_out_bound_i) return false;
  }

  std::array<Vector3<S>, 6> edges1;
  edges1[0] = vertices1[0] - vertices1[1];
  edges1[1] = vertices1[0] - vertices1[2];
  edges1[2] = vertices1[0] - vertices1[3];
  edges1[3] = vertices1[1] - vertices1[2];
  edges1[4] = vertices1[1] - vertices1[3];
  edges1[5] = vertices1[2] - vertices1[3];
  std::array<Vector3<S>, 6> edges2;
  edges2[0] = vertices2[0] - vertices2[1];
  edges2[1] = vertices2[0] - vertices2[2];
  edges2[2] = vertices2[0] - vertices2[3];
  edges2[3] = vertices2[1] - vertices2[2];
  edges2[4] = vertices2[1] - vertices2[3];
  edges2[5] = vertices2[2] - vertices2[3];

  std::array<Vector3<S>, 36> axis2test;
  for (int edge1 = 0, offset = 0; edge1 < 6; edge1++) {
    for (int edge2 = 0; edge2 < 6; edge2++) {
      axis2test[offset] = edges1[edge1].cross(edges2[edge2]);
      offset += 1;
    }
  }

  // Check separation along all axes
  for (const auto& axis : axis2test) {
    const bool is_separated =
        isSeparatedOnAxis<S, 4, 4>(vertices1, vertices2, axis);
    if (is_separated) return false;
  }

  return true;
}

template <typename S>
bool boxTerahedronIntersect(const Box<S>& s1, const Transform3<S>& tf1,
                            const Tetrahedron<S>& s2,
                            const Transform3<S>& tf2) {
  constexpr int n_box_vertices = 8;
  std::array<Vector3<S>, n_box_vertices> vertices1;
  {
    const auto a = s1.side[0] / 2;
    const auto b = s1.side[1] / 2;
    const auto c = s1.side[2] / 2;
    vertices1[0] = Vector3<S>(a, b, c);
    vertices1[1] = Vector3<S>(a, b, -c);
    vertices1[2] = Vector3<S>(a, -b, c);
    vertices1[3] = Vector3<S>(a, -b, -c);
    vertices1[4] = Vector3<S>(-a, b, c);
    vertices1[5] = Vector3<S>(-a, b, -c);
    vertices1[6] = Vector3<S>(-a, -b, c);
    vertices1[7] = Vector3<S>(-a, -b, -c);
  }
  const auto vertices2 =
      getTetrahedronVertices(s2, tf1.inverse(Eigen::Isometry) * tf2);
  std::array<Vector3<S>, kVertexOfTetrahedron> normals2;
  std::array<S, kVertexOfTetrahedron> height2;
  calcTetrahedronNormalHeight(vertices2, normals2, height2);

  // Check separation using tetrahedron faces and box vertices
  for (int i = 0; i < kVertexOfTetrahedron; i++) {
    const auto& normal_i = normals2[i];
    const auto& point_in_face_i = vertices2[(i + 1) % kVertexOfTetrahedron];

    const bool all_out_bound_i = checkPointDistanceToFaceAllOutsideBounds<S, 8>(
        vertices1, normal_i, point_in_face_i, height2[i], 0);
    if (all_out_bound_i) return false;
  }

  // Check separation using box faces and tetrahedron vertices
  for (int i = 0; i < 3; i++) {
    const auto face_i = s1.side[i] / 2;
    bool not_overlap =
        std::find_if_not(std::begin(vertices2), std::end(vertices2),
                         [&face_i, &i](const Vector3<S>& vertice) {
                           return vertice[i] > face_i;
                         }) == std::end(vertices2);
    if (not_overlap) return false;
    not_overlap = std::find_if_not(std::begin(vertices2), std::end(vertices2),
                                   [&face_i, &i](const Vector3<S>& vertice) {
                                     return vertice[i] < (-face_i);
                                   }) == std::end(vertices2);
    if (not_overlap) return false;
  }

  std::array<Vector3<S>, 6> edges_2;
  edges_2[0] = vertices2[0] - vertices2[1];
  edges_2[1] = vertices2[0] - vertices2[2];
  edges_2[2] = vertices2[0] - vertices2[3];
  edges_2[3] = vertices2[1] - vertices2[2];
  edges_2[4] = vertices2[1] - vertices2[3];
  edges_2[5] = vertices2[2] - vertices2[3];
  std::array<Vector3<S>, 3> edges_1{Vector3<S>(1, 0, 0), Vector3<S>(0, 1, 0),
                                    Vector3<S>(0, 0, 1)};

  std::array<Vector3<S>, 18> axis2test;
  int offset = 0;
  for (const auto& edge_1 : edges_1) {
    for (const auto& edge_2 : edges_2) {
      axis2test[offset] = edge_2.cross(edge_1);
      offset += 1;
    }
  }

  for (const auto& axis : axis2test) {
    const bool is_separated =
        isSeparatedOnAxis<S, 8, 4>(vertices1, vertices2, axis);
    if (is_separated) return false;
  }

  return true;
}

template <typename S>
bool triangleTerahedronIntersect(const TriangleP<S>& s1,
                                 const Transform3<S>& tf1,
                                 const Tetrahedron<S>& s2,
                                 const Transform3<S>& tf2) {
  std::array<Vector3<S>, kVertexOfTriangle> vertices1;
  vertices1[0] = tf1 * s1.a;
  vertices1[1] = tf1 * s1.b;
  vertices1[2] = tf1 * s1.c;
  const auto vertices2 = getTetrahedronVertices(s2, tf2);

  std::array<Vector3<S>, kVertexOfTetrahedron> normals2;
  std::array<S, kVertexOfTetrahedron> height2;
  calcTetrahedronNormalHeight(vertices2, normals2, height2);

  for (int i = 0; i < kVertexOfTetrahedron; i++) {
    const auto& normal_i = normals2[i];
    const auto& point_in_face_i = vertices2[(i + 1) % kVertexOfTetrahedron];

    const auto egde1 = vertices1[0] - point_in_face_i;
    const S distance1 = egde1.dot(normal_i);
    const auto egde2 = vertices1[1] - point_in_face_i;
    const S distance2 = egde2.dot(normal_i);
    const auto egde3 = vertices1[2] - point_in_face_i;
    const S distance3 = egde3.dot(normal_i);

    const bool in_positive_direction =
        distance1 > 0 && distance2 > 0 && distance3 > 0;
    if (in_positive_direction) return false;

    const bool in_negative_direction = distance1 < height2[i] &&
                                       distance2 < height2[i] &&
                                       distance3 < height2[i];
    if (in_negative_direction) return false;
  }

  std::array<Vector3<S>, 3> edges0;
  edges0[0] = vertices1[0] - vertices1[1];
  edges0[1] = vertices1[1] - vertices1[2];
  edges0[2] = vertices1[2] - vertices1[0];

  const auto tri_normal = edges0[0].cross(edges0[1]);
  {
    const bool all_out_bound = checkPointDistanceToFaceAllOutsideBounds<S, 4>(
        vertices2, tri_normal, vertices1[0], 0, 0);
    if (all_out_bound) return false;
  }

  std::array<Vector3<S>, kVertexOfTriangle> tri_edge_normals;
  std::array<S, 3> height1;
  for (int i = 0; i < kVertexOfTriangle; i++) {
    tri_edge_normals[i] = edges0[i].cross(tri_normal);
    const auto& v = vertices1[(i + 2) % kVertexOfTriangle];
    const auto edge = vertices1[i] - v;
    const S dot = edge.dot(tri_edge_normals[i]);

    if (dot > 0) {
      height1[i] = -dot;
    } else {
      tri_edge_normals[i] = -tri_edge_normals[i];
      height1[i] = dot;
    }
  }

  for (int i = 0; i < kVertexOfTriangle; i++) {
    const auto& normal_i = tri_edge_normals[i];
    const auto& point_in_face_i = vertices1[i];
    const bool all_out_bound_i = checkPointDistanceToFaceAllOutsideBounds<S, 4>(
        vertices2, normal_i, point_in_face_i, height1[i], 0);
    if (all_out_bound_i) return false;
  }

  std::array<Vector3<S>, 6> edges1;
  edges1[0] = vertices2[0] - vertices2[1];
  edges1[1] = vertices2[0] - vertices2[2];
  edges1[2] = vertices2[0] - vertices2[3];
  edges1[3] = vertices2[1] - vertices2[2];
  edges1[4] = vertices2[1] - vertices2[3];
  edges1[5] = vertices2[2] - vertices2[3];

  std::array<Vector3<S>, 18> axis2test;
  int offset = 0;
  for (const auto& edge0 : edges0) {
    for (const auto& edge1 : edges1) {
      axis2test[offset] = edge1.cross(edge0);
      offset += 1;
    }
  }

  // Check separation along all axes
  for (const auto& axis : axis2test) {
    const bool is_separated =
        isSeparatedOnAxis<S, 3, 4>(vertices1, vertices2, axis);
    if (is_separated) return false;
  }

  return true;
}
}  // namespace detail
}  // namespace fcl
