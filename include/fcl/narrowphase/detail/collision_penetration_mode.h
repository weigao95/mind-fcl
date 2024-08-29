//
// Created by mech-mind_gw on 3/16/2023.
//

#pragma once

#include "fcl/common/types.h"

namespace fcl {
namespace detail {

/// fcl::collide would perform collision detection, during which the
/// PENETRATION distance can be optionally computed for intersecting object
/// pairs. This option specified the algorithm & result for the penetration
/// distance computation.
enum class CollisionPenetrationType {
  /// The algorithm will only record pairs of geometry that intersect
  /// with each other. Their contact normal/center/penetration will NOT
  /// be computed.
  Disabled,

  /// The default option which compute minimum penetration distance. In
  /// other words, for the given interesting object pair the algorithm would
  /// find a direction such that the penetration distance is minimized.
  ///
  /// Note that penetration distance is accurate ONLY FOR CONVEX geometries,
  /// which include basic geometry and fcl::Convex. For bvh/octree/heightmap
  /// the penetration distance is actually a pairwise penetration distance
  /// among all its components (triangle for bvh, boxes for octree/heightmap)
  DefaultGJK_EPA,

  /// Compute the penetration direction along a given direction.
  /// The underline algorithm is a modified MPR.
  DirectedPenetration,

  /// Compute the minimum penetration distance with a hint direction, which
  /// might comes from spatial or temporal coherence. The algorithm starts
  /// from that direction and incremental update the direction as well as
  /// the penetration depth. This method can be very efficient especially
  /// a good direction hint is available.
  IncrementalMinimumPenetration
};

/// CollisionPenetrationType + unit Vector3 direction, as DirectedPenetration
/// and IncrementalMinimumPenetration request a prior direction.
template <typename S>
class CollisionPenetrationMode {
 public:
  // Default construct
  explicit CollisionPenetrationMode();
  explicit CollisionPenetrationMode(CollisionPenetrationType type);
  CollisionPenetrationMode(CollisionPenetrationType type,
                           Vector3<S> unit_direction);

  // Setting method
  void setAsDisabled();
  void setAsDefaultGJK_EPA();
  void setAsDirectedPenetration(Vector3<S> request_direction_unit);
  void setAsIncrementalMinimumPenetration(
      Vector3<S> direction_hint_unit);

  // Access method
  CollisionPenetrationType penetration_type() const { return type_; }
  const Vector3<S>& priorPenetrationDirection() const;

  // For commonly used types
  static const CollisionPenetrationMode<S>& Disabled();
  static const CollisionPenetrationMode<S>& MinimumPenetrationEPA();

 private:
  CollisionPenetrationType type_{CollisionPenetrationType::Disabled};
  Vector3<S> unit_direction_data_;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/collision_penetration_mode-inl.h"