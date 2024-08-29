//
// Created by wei on 22-6-5.
//

#pragma once

#include "vec_types.h"

namespace fcl {
namespace cvx_collide {

template <typename T>
struct MinDistanceToSimplex {
  bool witness_in_simplex{false};
  T min_distance_square;
  Vector3<T> witness;
};

template <typename T>
void pointToSegmentSquaredDistance(const Vector3<T>& p, const Vector3<T>& x0,
                                   const Vector3<T>& x1,
                                   bool& witness_in_segment,
                                   T& min_distance_square,
                                   Vector3<T>* witness) {
  Vector3<T> d = x1 - x0;
  Vector3<T> a = x0 - p;
  const T d_length_square = d.squaredNorm();
  T t = T(-1.0) * a.dot(d) / d_length_square;

  // Depends on cases
  if (t <= T(0)) {
    if (witness != nullptr) *witness = x0;
    witness_in_segment = false;
    min_distance_square = (x0 - p).squaredNorm();
  } else if (t >= T(1.0)) {
    if (witness != nullptr) *witness = x1;
    witness_in_segment = false;
    min_distance_square = (x1 - p).squaredNorm();
  } else {
    // reuse a?
    Vector3<T>& local_witness = a;
    local_witness = x0 + t * d;
    if (witness != nullptr) *witness = local_witness;
    witness_in_segment = true;
    min_distance_square = local_witness.squaredNorm();
  }
}

template <typename T>
void pointToSegmentSquaredDistance(const Vector3<T>& p, const Vector3<T>& x0,
                                   const Vector3<T>& x1,
                                   MinDistanceToSimplex<T>& output) {
  pointToSegmentSquaredDistance(p, x0, x1, output.witness_in_simplex,
                                output.min_distance_square, &output.witness);
}

template <typename T>
T pointToLineDistance(const Vector3<T>& p, const Vector3<T>& a,
                      const Vector3<T>& b) {
  const Vector3<T> a_to_b = b - a;
  const T ab_length = a_to_b.norm();
  if (ab_length <= T(0.0)) {
    return (p - a).norm();
  } else {
    const Vector3<T> plane_pab_normal = (p - a).cross(a_to_b);
    return plane_pab_normal.norm() / ab_length;
  }
}

template <typename T>
T pointToPlaneDistance(const Vector3<T>& p, const Vector3<T>& a,
                       const Vector3<T>& b, const Vector3<T>& c) {
  const Vector3<T> normal = (a - b).cross(b - c);
  const T normal_length = normal.norm();
  if (normal_length <= T(0.0)) {
    // Use the distance from p to LINE ab?
    return pointToLineDistance(p, a, b);
  }

  // No degenerate
  const Vector3<T> p_to_a = a - p;
  const Vector3<T> unit_normal = normal / normal_length;
  return std::abs(unit_normal.dot(p_to_a));
}

template <typename T>
void pointToTriangleSquaredDistance(const Vector3<T>& p, const Vector3<T>& a,
                                    const Vector3<T>& b, const Vector3<T>& c,
                                    bool& witness_in_triangle,
                                    T& min_distance_square,
                                    Vector3<T>* witness) {
  const Vector3<T> dl[] = {a - b, b - c, c - a};
  const Vector3<T> n = dl[0].cross(dl[1]);
  const T n_squared_norm = n.squaredNorm();
  const T area = std::sqrt(n_squared_norm);

  // Functor to handle sub-simplex
  auto process_distance_in_sub_simplex = [&]() -> void {
    std::array<MinDistanceToSimplex<T>, 3> distance_to_segments;
    pointToSegmentSquaredDistance<T>(p, a, b, distance_to_segments[0]);
    pointToSegmentSquaredDistance<T>(p, b, c, distance_to_segments[1]);
    pointToSegmentSquaredDistance<T>(p, c, a, distance_to_segments[2]);
    int min_distance_idx = 0;
    for (int i = 1; i < 3; i++) {
      if (distance_to_segments[i].min_distance_square <
          distance_to_segments[min_distance_idx].min_distance_square) {
        min_distance_idx = i;
      }
    }

    // Assign to output
    witness_in_triangle = false;
    min_distance_square =
        distance_to_segments[min_distance_idx].min_distance_square;
    if (witness) *witness = distance_to_segments[min_distance_idx].witness;
  };
  if (std::abs(area) <= T(1e-16)) {
    process_distance_in_sub_simplex();
    return;
  }

  // Vector from p to the PLANE-projection of p
  const T d = (a - p).dot(n);
  const Vector3<T> p_to_projected = n * (d / n_squared_norm);
  const Vector3<T> projected_p = p + p_to_projected;

  // Compute the parameter
  std::array<T, 3> parameterization;
  parameterization[0] = dl[1].cross(b - projected_p).norm() / area;
  parameterization[1] = dl[2].cross(c - projected_p).norm() / area;
  parameterization[2] = T(1.0) - parameterization[0] - parameterization[1];

  // Compute the parameterization in triangle
  witness_in_triangle = true;
  for (auto i = 0; i < 3; i++) {
    if (parameterization[i] < T(0.0) || parameterization[i] > T(1.0)) {
      witness_in_triangle = false;
      break;
    }
  }

  // Check with another method
  const T parameterization2_for_checking =
      dl[0].cross(a - projected_p).norm() / area;
  constexpr T difference_tol = T(1e-3);
  if (std::abs(parameterization2_for_checking - parameterization[2]) >
      difference_tol) {
    witness_in_triangle = false;
  }

  // In triangle
  if (witness_in_triangle) {
    min_distance_square = p_to_projected.squaredNorm();
    if (witness) {
      *witness = a * parameterization[0] + b * parameterization[1] +
                 c * parameterization[2];
    }

    // Done with in triangle case
    return;
  }

  // else, witness not in triangle
  process_distance_in_sub_simplex();
}

template <typename T>
void pointToTriangleSquaredDistance(const Vector3<T>& p, const Vector3<T>& a,
                                    const Vector3<T>& b, const Vector3<T>& c,
                                    MinDistanceToSimplex<T>& min_distance) {
  pointToTriangleSquaredDistance<T>(p, a, b, c, min_distance.witness_in_simplex,
                                    min_distance.min_distance_square,
                                    &min_distance.witness);
}

}  // namespace cvx_collide
}  // namespace fcl
