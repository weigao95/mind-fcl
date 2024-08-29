//
// Created by mech-mind_gw on 2/7/2022.
//

#pragma once

#include "minkowski_diff.h"
#include "vec_types.h"

namespace fcl {
namespace cvx_collide {

/// The Minkowski Portal Refinement Algorithm Described
/// in Game Programming Gem 7. This algorithm can be used
/// to computed binary intersection and DIRECTED penetration
template <typename T>
struct MPR {
  MPR(int max_iterations_in, T tolerance_in);
  ~MPR() = default;

  // The data for binary query
  enum class IntersectStatus { Intersect, Separated, Failed };
  struct IntersectData {
    Vector3<T> v0_interior;
    Vector3<T> v1;
    Vector3<T> v2;
    Vector3<T> v3;

    // The direction that find v1-v3
    Vector3<T> v1_dir_in_support;
    Vector3<T> v2_dir_in_support;
    Vector3<T> v3_dir_in_support;
  };

  // The evaluation interface with only binary flag output
  static IntersectStatus RunIntersect(const MinkowskiDiff<T>& shape,
                                      IntersectData& intersect_data,
                                      int max_iterations, T tolerance);
  IntersectStatus Intersect(const MinkowskiDiff<T>& shape,
                            IntersectData* intersect_data = nullptr) const;
  bool IsOriginEnclosedDebug(const MinkowskiDiff<T>& shape,
                             const IntersectData& intersect_data) const;

  // The data for penetration
  enum class DirectedPenetrationStatus {
    OK,
    Failed,
    FailedNoIntersect,
    FailedRefinementIterationLimit
  };
  struct DirectedPenetrationData {
    Vector3<T> v1;
    Vector3<T> v2;
    Vector3<T> v3;

    // The direction that find v1-v3
    Vector3<T> v1_dir_in_support;
    Vector3<T> v2_dir_in_support;
    Vector3<T> v3_dir_in_support;

    // The output
    Vector3<T> portal_normal;
    T distance_on_direction;
    Vector3<T> p0_in_shape0_frame;
    Vector3<T> p1_in_shape0_frame;
  };
  static DirectedPenetrationStatus RunDirectedPenetration(
      const MinkowskiDiff<T>& shape, const Vector3<T>& unit_direction,
      DirectedPenetrationData& penetration_data, int max_iterations,
      T tolerance);
  DirectedPenetrationStatus DirectedPenetration(
      const MinkowskiDiff<T>& shape, const Vector3<T>& unit_direction,
      DirectedPenetrationData& penetration_data) const;

  // The data for local refinement minimum penetration
  enum class IncrementalPenetrationStatus {
    OK,
    Failed,
    FailedNoIntersect,
    DirectedPenetrationFailed,
    IterationLimit
  };
  struct IncrementalMinimumPenetrationData {
    // The output
    T minimum_penetration;
    Vector3<T> penetration_direction;
    Vector3<T> p0_in_shape0_frame;
    Vector3<T> p1_in_shape0_frame;
  };
  static IncrementalPenetrationStatus RunIncrementalMinimumPenetrationDistance(
      const MinkowskiDiff<T>& shape, const Vector3<T>& init_direction,
      IncrementalMinimumPenetrationData& refine_data, int max_iteration,
      T tolerance, bool return_on_subroutine_converge = true);
  IncrementalPenetrationStatus IncrementalMinimumPenetrationDistance(
      const MinkowskiDiff<T>& shape, const Vector3<T>& init_direction,
      IncrementalMinimumPenetrationData& refine_data,
      bool return_on_subroutine_converge = true) const;

 private:
  // Termination
  const int max_iterations;
  const T tolerance;
  static constexpr T dot_eps_ratio = std::numeric_limits<T>::epsilon();

  // Helpers
  // Sometimes we need to normalize the direction (sphere, capsule, cylinder)
  // Sometimes we do not.
  static Vector3<T> computeShapeSupport(const MinkowskiDiff<T>& shape,
                                        Vector3<T>& direction);
  static inline T computeAbsNorm(const Vector3<T>& v);

  // Find the portal by finding v1/v2/v3 such that the origin ray from
  // v0 to O intersects with the triangle formed by v1/v2/v3.
  // v1/v2/v3 must be on the shape
  enum class FindPortalStatus { IterationLimit, DetectSeperated, PortalFound };
  static FindPortalStatus findPortal(
      const MinkowskiDiff<T>& shape, const Vector3<T>& v0, Vector3<T>& v1,
      Vector3<T>& v2, Vector3<T>& v3,
      std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support, int max_iterations);

  // Helpers for findPortal
  static void updatePortal(const Vector3<T>& v0, const Vector3<T>& v4,
                           const Vector3<T>& v4_dir_in_support, Vector3<T>& v1,
                           Vector3<T>& v2, Vector3<T>& v3,
                           std::array<Vector3<T>*, 3> v1_v2_v3_dir_in_support);

  static bool portalEncloseOrigin(const Vector3<T>& v0_interior,
                                  const Vector3<T>& v1, const Vector3<T>& v2,
                                  const Vector3<T>& v3);

  // Suppose there is a ray from origin O along the direction ray_direction
  // This ray intersects with the triangle v1v2v3, and the normal v123.dot(d)
  // cannot be negative. ray_direction/v123_normal is a unit vector.
  // Compute the distance from the origin to the intersecting point on v1v2v3.
  static void finalizeDirectedPenetrationResult(
      const MinkowskiDiff<T>& shape, const Vector3<T>& ray_direction,
      const Vector3<T>& v123_normal_arg,
      DirectedPenetrationData& penetration_data,
      const Vector3<T>* v4_optional = nullptr);
  static void finalizeIncrementalPenetrationResult(
      const MinkowskiDiff<T>& shape, const Vector3<T>& ray_direction,
      const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3,
      const Vector3<T>& v1_dir, const Vector3<T>& v2_dir,
      const Vector3<T>& v3_dir, const Vector3<T>& v123_normal_arg,
      T distance_lb, T distance_ub,
      IncrementalMinimumPenetrationData& penetration_data);

  // Private sub-routine for refinement
  enum class IncrementalMinDistanceSubroutineStatus {
    NewDirection,
    SubroutineConverge,
    Degenerated,
    Failed,
    FailedNoIntersect
  };
  static IncrementalMinDistanceSubroutineStatus
  incrementalMinimumDistanceExploreDirection(
      const MinkowskiDiff<T>& shape, const Vector3<T>& direction_to_explore,
      // Both input/output
      Vector3<T>& v1, Vector3<T>& v2, Vector3<T>& v3, Vector3<T>& v1_dir,
      Vector3<T>& v2_dir, Vector3<T>& v3_dir,
      // Output
      T& distance_on_direction_lb, T& distance_on_direction_ub,
      Vector3<T>& new_direction,
      // Parameters
      bool v123_valid, int max_iterations, T tolerance);
};

}  // namespace cvx_collide
}  // namespace fcl

// Header only implementation
#include "mpr.hpp"
#include "mpr_incremental_penetration.hpp"
