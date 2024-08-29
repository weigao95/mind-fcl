//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/common/types.h"
#include "fcl/geometry/collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/ccd_typedef.h"

namespace fcl {

template <typename S>
struct ContinuousCollisionContact {
  /// Collision geometry 1/2
  const CollisionGeometry<S>* o1{nullptr};
  const CollisionGeometry<S>* o2{nullptr};

  /// Same as usual contact
  static constexpr std::int64_t NONE = -1;
  std::int64_t b1{NONE};
  std::int64_t b2{NONE};

  /// For octree/heightmap collision
  /// If none of o1 and o2 is octree/heightmap, then ignore these results below
  /// If o1 and o2 are all octree/heightmap, then the name are clear
  /// If o1 OR o2 are octree/heightmap, then the bounding box is written to the
  /// corresponded objects (o1_bv or o2_bv).
  AABB<S> o1_bv{};
  AABB<S> o2_bv{};

  /// Time of collision interval in range of [0, 1], value
  /// might be invalid according to TocRequestType
  Interval<S> toc{};
};

}  // namespace fcl
