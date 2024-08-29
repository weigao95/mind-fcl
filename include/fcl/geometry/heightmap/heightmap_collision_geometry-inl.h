#pragma once

namespace fcl {

template <typename S>
HeightMapCollisionGeometry<S>::HeightMapCollisionGeometry(
    std::shared_ptr<const heightmap::LayeredHeightMap<S>> map)
    : height_map(std::move(map)) {}

template <typename S>
void HeightMapCollisionGeometry<S>::computeLocalAABB() {
  const S delta_x = height_map->half_range_x();
  const S delta_y = height_map->half_range_y();

  this->aabb_local = AABB<S>(
      Vector3<S>(-delta_x, -delta_y, 0),
      Vector3<S>(delta_x, delta_y, height_map->height_upper_bound_meter()));
  this->aabb_center = this->aabb_local.center();
  this->aabb_radius = (this->aabb_local.min_ - this->aabb_center).norm();
}

template <typename S>
OBJECT_TYPE HeightMapCollisionGeometry<S>::getObjectType() const {
  return OT_HEIGHTMAP;
}

template <typename S>
NODE_TYPE HeightMapCollisionGeometry<S>::getNodeType() const {
  return GEOM_HEIGHTMAP;
}

template <typename S>
const std::shared_ptr<const heightmap::LayeredHeightMap<S>>&
HeightMapCollisionGeometry<S>::raw_heightmap() const {
  return height_map;
}

}  // namespace fcl