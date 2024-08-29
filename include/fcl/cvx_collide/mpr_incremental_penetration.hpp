//
// Created by mech-mind_gw on 3/28/2022.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
void computeExploredDistanceLowerUpperBound(
    const Vector3<T>& ray_direction, const Vector3<T>& v1_arg,
    const Vector3<T>& v2_arg, const Vector3<T>& v3_arg,
    const Vector3<T>& v123_normal_arg, T o_to_v1_dot_v123_normal,
    T o_to_v4_dot_v123_normal, T& lower_bound, T& upper_bound) {
  // The input must be a UNIT direction
  const Vector3<T>& d = ray_direction;
  assert(std::abs(d.norm() - T(1.)) < T(1e-3));
  assert(std::abs(v123_normal_arg.norm() - T(1.)) < T(1e-3));
  const T v123_normal_dot_d = d.dot(v123_normal_arg);
  if (std::abs(v123_normal_dot_d) <= 0) {
    const T o_to_v1_dot_d = v1_arg.dot(d);
    const T o_to_v2_dot_d = v2_arg.dot(d);
    const T o_to_v3_dot_d = v3_arg.dot(d);
    lower_bound =
        std::max(o_to_v1_dot_d, std::max(o_to_v2_dot_d, o_to_v3_dot_d));
    upper_bound = lower_bound;
    return;
  }

  // The case we can make a division
  lower_bound = o_to_v1_dot_v123_normal / v123_normal_dot_d;
  upper_bound = o_to_v4_dot_v123_normal / v123_normal_dot_d;
  // assert(upper_bound >= lower_bound);
}

// The incremental minimum penetration query
template <typename T>
typename MPR<T>::IncrementalMinDistanceSubroutineStatus
MPR<T>::incrementalMinimumDistanceExploreDirection(
    const MinkowskiDiff<T>& shape, const Vector3<T>& direction_to_explore,
    // Both input/output
    Vector3<T>& v1, Vector3<T>& v2, Vector3<T>& v3,
    Vector3<T>& v1_dir_in_support, Vector3<T>& v2_dir_in_support,
    Vector3<T>& v3_dir_in_support,
    // Output
    T& distance_on_direction_lb, T& distance_on_direction_ub,
    Vector3<T>& new_direction,
    // Parameters
    bool v123_valid, int max_iterations, T tolerance) {
  using std::swap;
  const Vector3<T>& d = direction_to_explore;
  // v0_to_O is the direction. This v0 is actually a mock
  // It may NOT lie in the MinkowskiDiff shape
  const Vector3<T>& v0 = -d;

  // Compute the upper bound
  Vector3<T> init_support_dir = d;
  Vector3<T> init_support = computeShapeSupport(shape, init_support_dir);
  const T distance_ub_from_init_supporting = d.dot(init_support);

  // No intersect case
  if (distance_ub_from_init_supporting < 0) {
    return IncrementalMinDistanceSubroutineStatus::FailedNoIntersect;
  }

  // Do we need to re-init v123
  bool reinit_v123 = !v123_valid;

  // Init v1, v2 and v3
  if (reinit_v123) {
    v1_dir_in_support = d;
    v1 = init_support;

    v2_dir_in_support = v0.cross(v1);
    // o_to_v0 and o_to_v1 can be co-linear, check it
    // Note that v0 MUST have norm one
    if (computeAbsNorm(v2_dir_in_support) <= computeAbsNorm(v1) * tolerance) {
      // Refer to the note in RunIntersect
      distance_on_direction_lb = v1.dot(d);
      distance_on_direction_ub = distance_on_direction_lb;
      new_direction = d;
      return IncrementalMinDistanceSubroutineStatus::Degenerated;
    }

    v2 = computeShapeSupport(shape, v2_dir_in_support);
    if (v2_dir_in_support.dot(v2) < 0) {
      return IncrementalMinDistanceSubroutineStatus::FailedNoIntersect;
    }

    // it is better to form portal faces to be oriented "outside" origin
    // Here O must be an interior point in penetration query
    v3_dir_in_support = v1.cross(v2);
    if (v3_dir_in_support.dot(v0) > 0) {
      // swap v1/v2
      swap(v1, v2);
      swap(v1_dir_in_support, v2_dir_in_support);
      v3_dir_in_support *= -1;
    }

    v3 = computeShapeSupport(shape, v3_dir_in_support);
    if (v3_dir_in_support.dot(v3) < 0) {
      return IncrementalMinDistanceSubroutineStatus::FailedNoIntersect;
    }
  }

  // Scale the v0 to the max of v1/v2/v3
  Vector3<T> v0_scaled = v0;
  {
    const T v1_squared_norm = v1.squaredNorm();
    const T v2_squared_norm = v2.squaredNorm();
    const T v3_squared_norm = v3.squaredNorm();
    const T max_squared_norm =
        std::max(v1_squared_norm, std::max(v2_squared_norm, v3_squared_norm));
    v0_scaled *= std::sqrt(max_squared_norm);
  }

  // The loop to find the portal
  std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support{
      &v1_dir_in_support, &v2_dir_in_support, &v3_dir_in_support};
  auto find_portal_status = findPortal(shape, v0_scaled, v1, v2, v3,
                                       v1_v2_v3_dir_in_support, max_iterations);
  if (find_portal_status == FindPortalStatus::IterationLimit) {
    return IncrementalMinDistanceSubroutineStatus::Failed;
  } else if (find_portal_status == FindPortalStatus::DetectSeperated) {
    return IncrementalMinDistanceSubroutineStatus::FailedNoIntersect;
  } else {
    assert(find_portal_status == FindPortalStatus::PortalFound);
  }

  // Portal refinement
  int portal_refinement_iteration = 0;
  while (portal_refinement_iteration < max_iterations) {
    // Update iteration data
    portal_refinement_iteration += 1;

    // Compute the normal
    // The v123_normal must be oriented in the same side with O
    Vector3<T> v123_normal = (v2 - v1).cross(v3 - v1);
    if (v123_normal.dot(d) < 0) {
      swap(v2, v3);
      swap(v2_dir_in_support, v3_dir_in_support);
      v123_normal *= -1;
    }

    // A new point v4 on that direction
    // v123_normal would be normalized in computeShapeSupport
    const Vector3<T> v4 = computeShapeSupport(shape, v123_normal);
    assert(std::abs(v123_normal.norm() - T(1.)) < T(1e-3));
    const T distance_ub_on_new_supporting_plane = v4.dot(v123_normal);
    if (distance_ub_on_new_supporting_plane < 0) {
      return IncrementalMinDistanceSubroutineStatus::FailedNoIntersect;
    }

    // Compute some distance
    // const Vector3<T> v1_to_v4 = v4 - v1;
    // const T v1_v4_distance_on_n123 = v1_to_v4.dot(v123_normal);
    const T v1_dot_v123_normal = v1.dot(v123_normal);
    const T v4_dot_v123_normal = v4.dot(v123_normal);
    const T v1_v4_distance_on_n123 = v4_dot_v123_normal - v1_dot_v123_normal;
    computeExploredDistanceLowerUpperBound(
        d, v1, v2, v3, v123_normal, v1_dot_v123_normal, v4_dot_v123_normal,
        distance_on_direction_lb, distance_on_direction_ub);

    // Compute the shortcut
    if (distance_on_direction_lb >
        distance_ub_on_new_supporting_plane + tolerance) {
      new_direction = v123_normal;
      return IncrementalMinDistanceSubroutineStatus::NewDirection;
    }

    // Separation plane very close to the new (candidate) portal
    // Note that v123_normal is normalized
    // We must should stop the subroutine
    if (std::abs(v1_v4_distance_on_n123) <= tolerance) {
      new_direction = v123_normal;
      return IncrementalMinDistanceSubroutineStatus::SubroutineConverge;
    }

    // Update the portal
    updatePortal(v0_scaled, v4, v123_normal, v1, v2, v3,
                 v1_v2_v3_dir_in_support);
  }

  // Failure
  return IncrementalMinDistanceSubroutineStatus::Failed;
}

template <typename T>
typename MPR<T>::IncrementalPenetrationStatus
MPR<T>::RunIncrementalMinimumPenetrationDistance(
    const MinkowskiDiff<T>& shape, const Vector3<T>& init_direction,
    IncrementalMinimumPenetrationData& refine_data, int max_iteration,
    T tolerance, bool return_on_subroutine_converge) {
  // The local variable
  Vector3<T> d = init_direction;
  assert(std::abs(d.norm() - T(1.0)) <= T(1e-3));
  Vector3<T> v1, v2, v3;
  Vector3<T> v1_dir_in_support, v2_dir_in_support, v3_dir_in_support;

  // The main loop
  Vector3<T> prev_direction = init_direction;
  T prev_distance_lb{-1};
  T prev_distance_ub{-1};
  for (auto outer_iteration = 0; outer_iteration < max_iteration;
       outer_iteration++) {
    // Invoke subroutine
    Vector3<T> new_d;
    T distance_on_d_lb{0};
    T distance_on_d_ub{0};
    const bool portal_valid = (outer_iteration >= 1);
    auto subroutine_status = incrementalMinimumDistanceExploreDirection(
        shape, d, v1, v2, v3, v1_dir_in_support, v2_dir_in_support,
        v3_dir_in_support, distance_on_d_lb, distance_on_d_ub, new_d,
        portal_valid, max_iteration, tolerance);

    // Update the output
    switch (subroutine_status) {
      case IncrementalMinDistanceSubroutineStatus::NewDirection: {
        d = new_d;
        break;
      }
      case IncrementalMinDistanceSubroutineStatus::SubroutineConverge: {
        // Directly return without further checking.
        if (return_on_subroutine_converge) {
          finalizeIncrementalPenetrationResult(
              shape, d, v1, v2, v3, v1_dir_in_support, v2_dir_in_support,
              v3_dir_in_support, new_d, distance_on_d_lb, distance_on_d_ub,
              refine_data);
          return IncrementalPenetrationStatus::OK;
        }

        // Prev-data is valid, use that to check the convergence
        if (outer_iteration >= 1) {
          assert(prev_distance_lb >= 0);
          assert(prev_distance_ub >= 0);
          if (std::abs(distance_on_d_ub - prev_distance_ub) <= tolerance) {
            // Done as not much improvement
            finalizeIncrementalPenetrationResult(
                shape, d, v1, v2, v3, v1_dir_in_support, v2_dir_in_support,
                v3_dir_in_support, new_d, distance_on_d_lb, distance_on_d_ub,
                refine_data);
            return IncrementalPenetrationStatus::OK;
          }
        }

        // Check direction
        constexpr bool check_direction_termination = true;
        constexpr T direction_tolerance = T(1e-3);
        if (check_direction_termination &&
            ((new_d - prev_direction).norm() < direction_tolerance)) {
          finalizeIncrementalPenetrationResult(
              shape, d, v1, v2, v3, v1_dir_in_support, v2_dir_in_support,
              v3_dir_in_support, new_d, distance_on_d_lb, distance_on_d_ub,
              refine_data);
          return IncrementalPenetrationStatus::OK;
        }

        // Move to next, as the subroutine converge
        // The portal can be rather small here
        // portal_valid = false;
        d = new_d;
        break;
      }
      case IncrementalMinDistanceSubroutineStatus::Degenerated: {
        // Done on degenerate
        refine_data.minimum_penetration = distance_on_d_ub;
        refine_data.penetration_direction = new_d;
        refine_data.p0_in_shape0_frame = shape.support0(new_d);
        refine_data.p1_in_shape0_frame = shape.support1(-new_d);
        return IncrementalPenetrationStatus::OK;
      }
      case IncrementalMinDistanceSubroutineStatus::FailedNoIntersect: {
        return IncrementalPenetrationStatus::FailedNoIntersect;
      }
      default:
        return IncrementalPenetrationStatus::Failed;
    }

    // Move to next iteration, book keeping
    prev_direction = new_d;
    prev_distance_lb = distance_on_d_lb;
    prev_distance_ub = distance_on_d_ub;
  }

  // Too many iterations
  finalizeIncrementalPenetrationResult(shape, d, v1, v2, v3, v1_dir_in_support,
                                       v2_dir_in_support, v3_dir_in_support,
                                       prev_direction, prev_distance_lb,
                                       prev_distance_ub, refine_data);
  return IncrementalPenetrationStatus::IterationLimit;
}

template <typename T>
typename MPR<T>::IncrementalPenetrationStatus
MPR<T>::IncrementalMinimumPenetrationDistance(
    const MinkowskiDiff<T>& shape, const Vector3<T>& init_direction,
    IncrementalMinimumPenetrationData& refine_data,
    bool return_on_subroutine_converge) const {
  return RunIncrementalMinimumPenetrationDistance(
      shape, init_direction, refine_data, max_iterations, tolerance,
      return_on_subroutine_converge);
}

template <typename T>
void MPR<T>::finalizeIncrementalPenetrationResult(
    const MinkowskiDiff<T>& shape, const Vector3<T>& ray_direction,
    const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3,
    const Vector3<T>& v1_dir, const Vector3<T>& v2_dir,
    const Vector3<T>& v3_dir, const Vector3<T>& v123_normal_arg, T distance_lb,
    T distance_ub,
    typename cvx_collide::MPR<T>::IncrementalMinimumPenetrationData&
        penetration_data) {
  const Vector3<T>& d = ray_direction;
  penetration_data.minimum_penetration = distance_ub;
  penetration_data.penetration_direction = d;

  // Point in triangle distance is just lb
  (void)(v123_normal_arg);
  assert(std::abs(d.dot(v123_normal_arg)) > 0);
  Vector3<T> o_projected = distance_lb * d;
  const Vector3<T> v1_to_v2 = v2 - v1;
  const Vector3<T> v1_to_v3 = v3 - v1;
  const Vector3<T> s1s2s3_plane_normal = v1_to_v2.cross(v1_to_v3);
  const T area = s1s2s3_plane_normal.norm();
  const T s2_weight = (v1_to_v3.cross(v1 - o_projected)).norm() / area;
  const T s3_weight = (v1_to_v2.cross(v1 - o_projected)).norm() / area;
  const T s1_weight = T(1.0) - s2_weight - s3_weight;
  penetration_data.p0_in_shape0_frame = shape.support0(v1_dir) * s1_weight +
                                        shape.support0(v2_dir) * s2_weight +
                                        shape.support0(v3_dir) * s3_weight;
  penetration_data.p1_in_shape0_frame = shape.support1(-v1_dir) * s1_weight +
                                        shape.support1(-v2_dir) * s2_weight +
                                        shape.support1(-v3_dir) * s3_weight;
}

}  // namespace cvx_collide
}  // namespace fcl