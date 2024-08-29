#pragma once

#include "fcl/geometry/heightmap/layered_heightmap.h"

namespace fcl {

/// CollisionGeometry class for HeightMap for integration
/// in fcl::collide. Just a const shared_ptr to the actual
/// LayeredHeightMap.
template <typename S>
class HeightMapCollisionGeometry : public CollisionGeometry<S> {
 public:
  using HeightMapPtr = std::shared_ptr<const heightmap::LayeredHeightMap<S>>;
  explicit HeightMapCollisionGeometry(
      std::shared_ptr<const heightmap::LayeredHeightMap<S>> map);

  /// Compute the AABB<S> for the geometry in its local coordinate system
  void computeLocalAABB() override;

  /// Simple access
  OBJECT_TYPE getObjectType() const override;
  NODE_TYPE getNodeType() const override;
  const HeightMapPtr& raw_heightmap() const;

  /// Read-only, shared_ptr access to a height_map
 private:
  std::shared_ptr<const heightmap::LayeredHeightMap<S>> height_map;
};

// template specialized for float/double
using HeightMapCollisionGeometryf = HeightMapCollisionGeometry<float>;
using HeightMapCollisionGeometryd = HeightMapCollisionGeometry<double>;

}  // namespace fcl

#include "fcl/geometry/heightmap/heightmap_collision_geometry-inl.h"