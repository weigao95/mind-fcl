#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
MPR<T>::MPR(int max_iterations_in, T tolerance_in)
    : max_iterations(max_iterations_in), tolerance(tolerance_in) {}

template <typename T>
Vector3<T> MPR<T>::computeShapeSupport(const MinkowskiDiff<T>& shape,
                                       Vector3<T>& direction) {
  direction.normalize();
  return shape.support(direction);
}

template <typename T>
T MPR<T>::computeAbsNorm(const Vector3<T>& v) {
  return std::abs(v[0]) + std::abs(v[1]) + std::abs(v[2]);
}

template <typename T>
typename MPR<T>::IntersectStatus MPR<T>::Intersect(
    const MinkowskiDiff<T>& shape, IntersectData* intersect_data) const {
  if (intersect_data != nullptr) {
    return RunIntersect(shape, *intersect_data, max_iterations, tolerance);
  } else {
    IntersectData data;
    return RunIntersect(shape, data, max_iterations, tolerance);
  }
}

template <typename T>
typename MPR<T>::IntersectStatus MPR<T>::RunIntersect(
    const MinkowskiDiff<T>& shape, IntersectData& intersect_data,
    int max_iterations, T tolerance) {
  // Gather the data
  using std::swap;
  Vector3<T>& v0_interior = intersect_data.v0_interior;
  Vector3<T>& v1 = intersect_data.v1;
  Vector3<T>& v2 = intersect_data.v2;
  Vector3<T>& v3 = intersect_data.v3;
  Vector3<T>& v1_dir_in_support = intersect_data.v1_dir_in_support;
  Vector3<T>& v2_dir_in_support = intersect_data.v2_dir_in_support;
  Vector3<T>& v3_dir_in_support = intersect_data.v3_dir_in_support;

  // Init the data as NaN
  v1.setConstant(std::numeric_limits<T>::quiet_NaN());
  v2.setConstant(std::numeric_limits<T>::quiet_NaN());
  v3.setConstant(std::numeric_limits<T>::quiet_NaN());
  assert(v1.array().isNaN().all());
  assert(v2.array().isNaN().all());
  assert(v3.array().isNaN().all());

  // Compute the interior point
  v0_interior = shape.interior();

  // Distance is smaller than a threshold, then must intersect
  if (v0_interior.squaredNorm() <= tolerance * tolerance) {
    return IntersectStatus::Intersect;
  }

  // Init v1, v2 and v3
  v1_dir_in_support = -v0_interior;
  v1 = computeShapeSupport(shape, v1_dir_in_support);

  // No intersect case 1
  if (v1_dir_in_support.dot(v1) < 0) {
    return IntersectStatus::Separated;
  }

  v2_dir_in_support = v0_interior.cross(v1);
  // o_to_v0 and o_to_v1 can be co-linear, check it
  // Note that v0 MUST have norm one
  // This equation might be written as
  //   cross(o_to_v0.normalized(), o_to_v1.normalized()).norm() <= tolerance
  // which can be further expanded as
  //   cross(o_to_v0, o_to_v1).norm() <= tolerance * v0.norm() * v1.norm()
  // However, we do not want to compute the L2 norm, thus replace it with
  // abs norm, which is much easier to compute.
  if (computeAbsNorm(v2_dir_in_support) <=
      computeAbsNorm(v0_interior) * computeAbsNorm(v1) * tolerance) {
    // o_to_v0 and o_to_v1 can be co-linear, from the condition above
    // v1_dir_in_support.dot(v1) = - v0_interior.dot(v1) < 0 is False
    // which implies
    // - v0_interior.dot(v1) > 0  --> v0_interior.dot(v1) < 0
    // As v0 = o_to_v0, v1 = o_to_v1, we have
    // o_to_v0.dot(o_to_v1) < 0, o is in the middle of v0/v1
    // As v0 is an interior point, v1 is a boundary point
    // We conclude O must be within the shape
    return IntersectStatus::Intersect;
  }

  v2 = computeShapeSupport(shape, v2_dir_in_support);
  if (v2_dir_in_support.dot(v2) < 0) {
    return IntersectStatus::Separated;
  }

  v3_dir_in_support = (v1 - v0_interior).cross(v2 - v0_interior);
  // it is better to form portal faces to be oriented "outside" origin
  if (v3_dir_in_support.dot(v0_interior) > 0) {
    // swap v1/v2
    swap(v1, v2);
    swap(v1_dir_in_support, v2_dir_in_support);
    v3_dir_in_support *= -1;
  }

  v3 = computeShapeSupport(shape, v3_dir_in_support);
  if (v3_dir_in_support.dot(v3) < 0) {
    return IntersectStatus::Separated;
  }

  // The loop to find the portal
  std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support{
      &v1_dir_in_support, &v2_dir_in_support, &v3_dir_in_support};
  auto find_portal_status = findPortal(shape, v0_interior, v1, v2, v3,
                                       v1_v2_v3_dir_in_support, max_iterations);
  if (find_portal_status == FindPortalStatus::IterationLimit) {
    return IntersectStatus::Failed;
  } else if (find_portal_status == FindPortalStatus::DetectSeperated) {
    return IntersectStatus::Separated;
  } else {
    assert(find_portal_status == FindPortalStatus::PortalFound);
  }

  // Portal refinement
  int portal_refinement_iteration = 0;
  while (portal_refinement_iteration < max_iterations) {
    // Update iteration data
    portal_refinement_iteration += 1;
    const Vector3<T>& o_to_v1 = v1;  // - o
    const Vector3<T>& o_to_v0 = v0_interior;

    // Compute the normal
    // The v123_normal must be oriented in the same side with O
    Vector3<T> v123_normal = (v2 - v1).cross(v3 - v1);
    // v123_normal should be in the same direction as v0_to_o
    // or the center ray
    if (v123_normal.dot(o_to_v0) > 0) {
      swap(v2, v3);
      swap(v2_dir_in_support, v3_dir_in_support);
      v123_normal *= -1;
    }

    // Check intersection
    {
      // Old impl
      // const bool v123_seperated_v0_O =
      //     is_sign_matched(o_to_v1.dot(v123_normal),
      //     v1_to_v0.dot(v123_normal));
      const T v123n_dot_o_to_v1 = o_to_v1.dot(v123_normal);
      // const Vector3<T> v1_to_v0 = v0_interior - v1;
      // const T v123n_dot_v1_to_v0 = v1_to_v0.dot(v123_normal);
      // It can be proved that dot(v123n, v0_to_v1) must be positive in the
      // convention mentioned above.
      // assert(v123n_dot_v1_to_v0 <= 0);
      const bool v123_seperated_v0_O = v123n_dot_o_to_v1 < 0;
      if (!v123_seperated_v0_O) {
        return IntersectStatus::Intersect;
      }
    }

    // A new point v4 on that direction
    const Vector3<T> v4 = computeShapeSupport(shape, v123_normal);
    if (v4.dot(v123_normal) < 0) {
      return IntersectStatus::Separated;
    }

    // Separation plane very close to the new (candidate) portal
    // Note that v123_normal can be un-normalized, thus its length
    // must be considered
    const Vector3<T> v1_to_v4 = v4 - v1;
    if (std::abs(v1_to_v4.dot(v123_normal)) <
        tolerance * computeAbsNorm(v123_normal)) {
      return IntersectStatus::Separated;
    }

    // Update the portal
    updatePortal(v0_interior, v4, v123_normal, v1, v2, v3,
                 v1_v2_v3_dir_in_support);
  }

  // Default case
  return IntersectStatus::Failed;
}

template <typename T>
typename MPR<T>::FindPortalStatus MPR<T>::findPortal(
    const MinkowskiDiff<T>& shape, const Vector3<T>& v0, Vector3<T>& v1,
    Vector3<T>& v2, Vector3<T>& v3,
    std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support, int max_iterations) {
  // Check input
  assert(!v1.array().isNaN().any());
  assert(!v2.array().isNaN().any());
  assert(!v3.array().isNaN().any());

  // The output buffer
  using std::swap;
  const Vector3<T>& o_to_v0 = v0;
  const T v0_abs_norm = computeAbsNorm(v0);
  Vector3<T>* v1_dir_in_support = v1_v2_v3_dir_in_support[0];
  Vector3<T>* v2_dir_in_support = v1_v2_v3_dir_in_support[1];
  Vector3<T>* v3_dir_in_support = v1_v2_v3_dir_in_support[2];
  int find_candidate_portal_iteration = 0;

  // The actual loop
  while (true) {
    if (find_candidate_portal_iteration >= max_iterations) {
      return FindPortalStatus::IterationLimit;
    }

    // Update iteration data
    find_candidate_portal_iteration += 1;
    Vector3<T> v0_to_v1 = v1 - v0;
    Vector3<T> v0_to_v2 = v2 - v0;
    Vector3<T> v0_to_v3 = v3 - v0;

    // Update the corresponded vertex
    // These normal are not oriented
    Vector3<T> v031_normal = v0_to_v3.cross(v0_to_v1);
    Vector3<T> v012_normal = v0_to_v1.cross(v0_to_v2);
    const T signed_volume = v0_to_v2.dot(v031_normal);
    // Orient it
    if (signed_volume < 0) {
      swap(v2, v3);
      swap(*v2_dir_in_support, *v3_dir_in_support);
      swap(v0_to_v2, v0_to_v3);

      // Something tricky here, note the changing of vectors
      swap(v012_normal, v031_normal);
      v031_normal *= -1;
      v012_normal *= -1;
      // signed_volume *= -1;
    }

    // Old impl
    // const bool v031_seperated_v2_and_o =
    //     is_sign_matched(o_to_v0.dot(v031_normal), v0_to_v2.dot(v031_normal));
    // New impl after orienting the v0123:
    assert(v0_to_v2.dot(v031_normal) >= 0);
    const bool v031_seperated_v2_and_o =
        (o_to_v0.dot(v031_normal) >
         dot_eps_ratio * v0_abs_norm * computeAbsNorm(v031_normal));
    if (v031_seperated_v2_and_o) {
      // Orient the normal towards O
      assert(o_to_v0.dot(v031_normal) > 0);
      Vector3<T>& search_v2_dir = v031_normal;
      search_v2_dir *= -1;

      // Find a new v2 in that direction
      v2 = computeShapeSupport(shape, search_v2_dir);
      if (v2_dir_in_support != nullptr) {
        *v2_dir_in_support = search_v2_dir;
      }

      // Miss detection
      if (v2.dot(search_v2_dir) < 0) {
        return FindPortalStatus::DetectSeperated;
      }

      continue;
    }

    // case 012
    // Old impl:
    // const bool v012_seperated_v3_and_o =
    //     is_sign_matched(o_to_v0.dot(v012_normal), v0_to_v3.dot(v012_normal));
    // New impl after orienting the v0123:
    assert(v0_to_v3.dot(v012_normal) >= 0);
    const bool v012_seperated_v3_and_o =
        (o_to_v0.dot(v012_normal) >
         dot_eps_ratio * v0_abs_norm * computeAbsNorm(v012_normal));
    if (v012_seperated_v3_and_o) {
      // Orient the normal towards O
      assert(o_to_v0.dot(v012_normal) > 0);
      Vector3<T>& search_v3_dir = v012_normal;
      search_v3_dir *= -1;

      // Find a new v3 in that direction
      v3 = computeShapeSupport(shape, search_v3_dir);
      if (v3_dir_in_support != nullptr) {
        *v3_dir_in_support = search_v3_dir;
      }

      // Miss detection
      if (v3.dot(search_v3_dir) < 0) {
        return FindPortalStatus::DetectSeperated;
      }

      // Loop again
      continue;
    }

    // case 023
    Vector3<T> v023_normal = v0_to_v2.cross(v0_to_v3);
    // Old impl
    // const bool v023_seperated_v1_and_o =
    //     is_sign_matched(o_to_v0.dot(v023_normal), v0_to_v1.dot(v023_normal));
    // New impl after orienting the v0123:
    assert(v0_to_v1.dot(v023_normal) >= 0);
    const bool v023_seperated_v1_and_o =
        (o_to_v0.dot(v023_normal) >
         dot_eps_ratio * v0_abs_norm * computeAbsNorm(v023_normal));
    if (v023_seperated_v1_and_o) {
      // Orient the normal towards O
      assert(o_to_v0.dot(v023_normal) > 0);
      Vector3<T>& search_v1_dir = v023_normal;
      search_v1_dir *= -1;

      // Find a new v2 in that direction
      v1 = computeShapeSupport(shape, search_v1_dir);
      if (v1_dir_in_support != nullptr) {
        *v1_dir_in_support = search_v1_dir;
      }

      // Miss detection
      if (v1.dot(search_v1_dir) < 0) {
        return FindPortalStatus::DetectSeperated;
      }

      // Loop again
      continue;
    }

    // No separation, we're done
    return FindPortalStatus::PortalFound;
  }
}

template <typename T>
bool MPR<T>::portalEncloseOrigin(const Vector3<T>& v0_interior,
                                 const Vector3<T>& v1, const Vector3<T>& v2,
                                 const Vector3<T>& v3) {
  // This is a debug method, thus we might place more attention on correctness
  using std::swap;
  const Vector3<T>& o_to_v0 = v0_interior;

  // No degeneration in this method
  assert(!v0_interior.array().isNaN().any());
  assert(!v1.array().isNaN().any());
  assert(!v2.array().isNaN().any());
  assert(!v3.array().isNaN().any());

  // Check the three faces contains O
  {
    Vector3<T> v0_to_v1 = v1 - v0_interior;
    Vector3<T> v0_to_v2 = v2 - v0_interior;
    Vector3<T> v0_to_v3 = v3 - v0_interior;

    // Update the corresponded vertex
    // These normal are not oriented
    Vector3<T> v031_normal = v0_to_v3.cross(v0_to_v1);
    Vector3<T> v012_normal = v0_to_v1.cross(v0_to_v2);
    const T signed_volume = v0_to_v2.dot(v031_normal);
    // Orient it
    if (signed_volume < 0) {
      swap(v0_to_v2, v0_to_v3);

      // Something tricky here, note the changing of vectors
      swap(v012_normal, v031_normal);
      v031_normal *= -1;
      v012_normal *= -1;
      // signed_volume *= -1;
    }

    // After orientation
    // Note that v2 and v3 might be swapped
    Vector3<T> v023_normal = v0_to_v2.cross(v0_to_v3);
    assert(v0_to_v2.dot(v031_normal) >= 0);
    assert(v0_to_v3.dot(v012_normal) >= 0);
    assert(v0_to_v1.dot(v023_normal) >= 0);

    // Compute separation
    const T v0_abs_norm = computeAbsNorm(o_to_v0);
    const bool v013_seperated_v2_and_o =
        o_to_v0.dot(v031_normal) >
        dot_eps_ratio * v0_abs_norm * computeAbsNorm(v031_normal);
    const bool v012_seperated_v3_and_o =
        o_to_v0.dot(v012_normal) >
        dot_eps_ratio * v0_abs_norm * computeAbsNorm(v012_normal);
    const bool v023_seperated_v1_and_o =
        o_to_v0.dot(v023_normal) >
        dot_eps_ratio * v0_abs_norm * computeAbsNorm(v023_normal);
    if (v013_seperated_v2_and_o || v012_seperated_v3_and_o ||
        v023_seperated_v1_and_o)
      return false;
  }

  // Check the last face without O, the face v123
  const Vector3<T>& o_to_v1 = v1;  // - o
  const Vector3<T> v1_to_v0 = v0_interior - v1;
  const Vector3<T> v1_to_v2 = v2 - v1;
  const Vector3<T> v1_to_v3 = v3 - v1;
  const Vector3<T> v123_normal = v1_to_v2.cross(v1_to_v3);
  const T v123n_dot_o_to_v1 = o_to_v1.dot(v123_normal);
  const T v123n_dot_v1_to_v0 = v1_to_v0.dot(v123_normal);
  const bool v123_seperated_v0_O =
      (v123n_dot_o_to_v1 > 0 && v123n_dot_v1_to_v0 > 0) ||
      (v123n_dot_o_to_v1 < 0 && v123n_dot_v1_to_v0 < 0);
  if (v123_seperated_v0_O) return false;

  // All faces have been checked
  return true;
}

template <typename T>
bool MPR<T>::IsOriginEnclosedDebug(const MinkowskiDiff<T>& shape,
                                   const IntersectData& intersect_data) const {
  // Degeneration case 0: one point
  assert(!intersect_data.v0_interior.array().isNaN().any());
  auto v0_local = shape.interior();
  if (intersect_data.v1.array().isNaN().any()) {
    return (v0_local.squaredNorm() <= tolerance * tolerance);
  }

  // Degeneration case 1: one segment
  assert(!intersect_data.v1.array().isNaN().any());
  auto v1_local = shape.support(intersect_data.v1_dir_in_support);
  if (intersect_data.v2.array().isNaN().any()) {
    return (computeAbsNorm(intersect_data.v2_dir_in_support) <=
            computeAbsNorm(v0_local) * computeAbsNorm(v1_local) * tolerance);
  }

  // Directly report not enclosed upon v3 NaN
  if (intersect_data.v3.array().isNaN().any()) {
    return false;
  }

  // Other points
  auto v2_local = shape.support(intersect_data.v2_dir_in_support);
  auto v3_local = shape.support(intersect_data.v3_dir_in_support);
  return portalEncloseOrigin(v0_local, v1_local, v2_local, v3_local);
}

template <typename T>
void MPR<T>::updatePortal(const Vector3<T>& v0, const Vector3<T>& v4,
                          const Vector3<T>& v123_normal, Vector3<T>& v1,
                          Vector3<T>& v2, Vector3<T>& v3,
                          std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support) {
  // Check input
  assert(!v1.array().isNaN().any());
  assert(!v2.array().isNaN().any());
  assert(!v3.array().isNaN().any());
  assert(!v4.array().isNaN().any());

  // The output buffer
  Vector3<T>* v1_dir_in_support = v1_v2_v3_dir_in_support[0];
  Vector3<T>* v2_dir_in_support = v1_v2_v3_dir_in_support[1];
  Vector3<T>* v3_dir_in_support = v1_v2_v3_dir_in_support[2];

  // v4 must appear in the next portal
  // select two in v1, v2 and v3
  // First do a separation in with plane v0_v4_o
  const Vector3<T>& o_to_v4 = v4;
  const Vector3<T>& o_to_v0 = v0;
  const Vector3<T>& o_to_v1 = v1;
  const Vector3<T>& o_to_v2 = v2;
  const Vector3<T>& o_to_v3 = v3;

  const Vector3<T> v0_v4_o_normal = o_to_v4.cross(o_to_v0);
  T dot = o_to_v1.dot(v0_v4_o_normal);
  if (dot > 0) {
    dot = o_to_v2.dot(v0_v4_o_normal);
    if (dot > 0) {
      // Discard v1
      v1 = v4;
      if (v1_dir_in_support != nullptr) {
        *v1_dir_in_support = v123_normal;
      }
    } else {
      // Discard v3
      v3 = v4;
      if (v3_dir_in_support != nullptr) {
        *v3_dir_in_support = v123_normal;
      }
    }
  } else {
    dot = o_to_v3.dot(v0_v4_o_normal);
    if (dot > 0) {
      // Discard v2
      v2 = v4;
      if (v2_dir_in_support != nullptr) {
        *v2_dir_in_support = v123_normal;
      }
    } else {
      // Discard v1
      v1 = v4;
      if (v1_dir_in_support != nullptr) {
        *v1_dir_in_support = v123_normal;
      }
    }
  }
}

// The directed penetration query
template <typename T>
typename MPR<T>::DirectedPenetrationStatus MPR<T>::RunDirectedPenetration(
    const MinkowskiDiff<T>& shape, const Vector3<T>& unit_direction,
    DirectedPenetrationData& penetration_data, int max_iterations,
    T tolerance) {
  using std::swap;
  const Vector3<T>& d = unit_direction;
  // v0_to_O is the direction. This v0 is actually a mock
  // It may NOT lie in the MinkowskiDiff shape
  const Vector3<T>& v0 = -d;

  // Gather the data
  Vector3<T>& v1 = penetration_data.v1;
  Vector3<T>& v2 = penetration_data.v2;
  Vector3<T>& v3 = penetration_data.v3;
  Vector3<T>& v1_dir_in_support = penetration_data.v1_dir_in_support;
  Vector3<T>& v2_dir_in_support = penetration_data.v2_dir_in_support;
  Vector3<T>& v3_dir_in_support = penetration_data.v3_dir_in_support;

  // Init v1, v2 and v3
  v1_dir_in_support = d;
  v1 = computeShapeSupport(shape, v1_dir_in_support);

  // No intersect case 1
  if (v1_dir_in_support.dot(v1) < 0) {
    return DirectedPenetrationStatus::FailedNoIntersect;
  }

  v2_dir_in_support = v0.cross(v1);
  // o_to_v0 and o_to_v1 can be co-linear, check it
  // Note that v0 MUST have norm one
  if (computeAbsNorm(v2_dir_in_support) <= computeAbsNorm(v1) * tolerance) {
    // Refer to the note in RunIntersect
    penetration_data.distance_on_direction = v1.dot(d);
    penetration_data.portal_normal = v1_dir_in_support;
    penetration_data.p0_in_shape0_frame = shape.support0(v1_dir_in_support);
    penetration_data.p1_in_shape0_frame = shape.support1(-v1_dir_in_support);
    return DirectedPenetrationStatus::OK;
  }

  v2 = computeShapeSupport(shape, v2_dir_in_support);
  if (v2_dir_in_support.dot(v2) < 0) {
    return DirectedPenetrationStatus::FailedNoIntersect;
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
    return DirectedPenetrationStatus::FailedNoIntersect;
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
    return DirectedPenetrationStatus::Failed;
  } else if (find_portal_status == FindPortalStatus::DetectSeperated) {
    return DirectedPenetrationStatus::FailedNoIntersect;
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
    const Vector3<T> v4 = computeShapeSupport(shape, v123_normal);
    if (v4.dot(v123_normal) < 0) {
      return DirectedPenetrationStatus::FailedNoIntersect;
    }

    // Separation plane very close to the new (candidate) portal
    // Note that v123_normal can be un-normalized, thus its length
    // must be considered
    const Vector3<T> v1_to_v4 = v4 - v1;
    if (std::abs(v1_to_v4.dot(v123_normal)) <
        tolerance * computeAbsNorm(v123_normal)) {
      finalizeDirectedPenetrationResult(shape, d, v123_normal, penetration_data,
                                        &v4);
      return DirectedPenetrationStatus::OK;
    }

    // Update the portal
    updatePortal(v0_scaled, v4, v123_normal, v1, v2, v3,
                 v1_v2_v3_dir_in_support);
  }

  // Iteration limit
  return DirectedPenetrationStatus::FailedRefinementIterationLimit;
}

template <typename T>
typename MPR<T>::DirectedPenetrationStatus MPR<T>::DirectedPenetration(
    const MinkowskiDiff<T>& shape, const Vector3<T>& unit_direction,
    DirectedPenetrationData& penetration_data) const {
  return RunDirectedPenetration(shape, unit_direction, penetration_data,
                                max_iterations, tolerance);
}

template <typename T>
void MPR<T>::finalizeDirectedPenetrationResult(
    const MinkowskiDiff<T>& shape, const Vector3<T>& ray_direction,
    const Vector3<T>& v123_normal_arg,
    DirectedPenetrationData& penetration_data, const Vector3<T>* v4_optional) {
  // First assign the normal
  penetration_data.portal_normal = v123_normal_arg;

  // Other terms require computation
  const Vector3<T>& d = ray_direction;
  const Vector3<T>& v1 = penetration_data.v1;
  const Vector3<T>& v2 = penetration_data.v2;
  const Vector3<T>& v3 = penetration_data.v3;
  const T v123_normal_dot_d = d.dot(v123_normal_arg);

  // Very unlikely case that can not divide the dot
  if (std::abs(v123_normal_dot_d) <= 0) {
    std::array<T, 3> d_dot_v123;
    std::array<Vector3<T>*, 3> d_for_v123;
    d_dot_v123[0] = v1.dot(d);
    d_dot_v123[1] = v2.dot(d);
    d_dot_v123[2] = v3.dot(d);
    d_for_v123[0] = &(penetration_data.v1_dir_in_support);
    d_for_v123[1] = &(penetration_data.v2_dir_in_support);
    d_for_v123[2] = &(penetration_data.v3_dir_in_support);
    T max_distance = -std::numeric_limits<T>::infinity();
    int max_distance_i = 0;
    for (int i = 0; i < 3; i++) {
      if (d_dot_v123[i] > max_distance) {
        max_distance = d_dot_v123[i];
        max_distance_i = i;
      }
    }

    // Assign the distance
    penetration_data.distance_on_direction = max_distance;

    // Assign the point
    const Vector3<T>& d_max_distance = *d_for_v123[max_distance_i];
    penetration_data.p0_in_shape0_frame = shape.support0(d_max_distance);
    penetration_data.p1_in_shape0_frame = shape.support1(-d_max_distance);
    return;
  }

  // The case we can make a division
  const T o_to_v_dot_v123_normal = v1.dot(v123_normal_arg);
  assert(o_to_v_dot_v123_normal >= 0);
  const T distance_to_v123 = o_to_v_dot_v123_normal / v123_normal_dot_d;
  if (v4_optional == nullptr) {
    penetration_data.distance_on_direction = distance_to_v123;
  } else {
    penetration_data.distance_on_direction =
        v4_optional->dot(v123_normal_arg) / v123_normal_dot_d;
  }

  // Point in triangle
  Vector3<T> o_projected = distance_to_v123 * d;
  const Vector3<T> v1_to_v2 = v2 - v1;
  const Vector3<T> v1_to_v3 = v3 - v1;
  const Vector3<T> s1s2s3_plane_normal = v1_to_v2.cross(v1_to_v3);
  const T area = s1s2s3_plane_normal.norm();
  const T s2_weight = (v1_to_v3.cross(v1 - o_projected)).norm() / area;
  const T s3_weight = (v1_to_v2.cross(v1 - o_projected)).norm() / area;
  const T s1_weight = T(1.0) - s2_weight - s3_weight;
  penetration_data.p0_in_shape0_frame =
      shape.support0(penetration_data.v1_dir_in_support) * s1_weight +
      shape.support0(penetration_data.v2_dir_in_support) * s2_weight +
      shape.support0(penetration_data.v3_dir_in_support) * s3_weight;
  penetration_data.p1_in_shape0_frame =
      shape.support1(-penetration_data.v1_dir_in_support) * s1_weight +
      shape.support1(-penetration_data.v2_dir_in_support) * s2_weight +
      shape.support1(-penetration_data.v3_dir_in_support) * s3_weight;
}

}  // namespace cvx_collide
}  // namespace fcl
