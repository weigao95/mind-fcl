//
// Created by mech-mind_gw on 2/22/2023.
//

#pragma once

#include "fcl/geometry/heightmap/flat_heightmap.h"

namespace fcl {
namespace heightmap {

/// LayeredHeightMap represents the same geometry as FlatHeightMap,
/// but it uses a stack of FlatHeightMaps for efficient collision
/// detection, especially with bvh/octree/itself.
/// The bottom layer is the one with the largest shape that
/// represents the finest geometry. For each layer except the bottom,
/// its one pixel represents a 2x2 grid of its lower layer.
///
/// Same coordinate and pixel convention is used in each layer.
/// However, pixel in each layer has difference RANGE.
template <typename S>
class LayeredHeightMap {
 private:
  /// The stack of layer. layers_[0] is top (the most coarse one)
  /// while layers_.back() is the bottom.
  static constexpr std::size_t kTopLayerIndex = 0;
  std::vector<FlatHeightMap<S>> layers_;
  FlatHeightMap<S>& bottom_mutable() { return layers_.back(); }

 public:
  /// The same construction parameter as FlatHeightMap,
  /// However, half_map_shape must be 2^n.
  LayeredHeightMap(S bottom_resolution_x, S bottom_resolution_y,
                   uint16_t bottom_half_map_shape_x,
                   uint16_t bottom_half_map_shape_y);
  LayeredHeightMap(S bottom_resolution, uint16_t bottom_half_map_shape);
  // Default copy/move/assign

  // The access of layers
  std::size_t n_layers() const { return layers_.size(); }
  const FlatHeightMap<S>& top() const { return layers_.front(); }
  const FlatHeightMap<S>& bottom() const { return layers_.back(); }
  const std::vector<FlatHeightMap<S>>& layers() const { return layers_; }
  bool is_bottom_layer(std::size_t layer_idx) const;
  // These information are the same for each layer
  uint16_t height_upper_bound_mm() const;
  S height_upper_bound_meter() const;
  S half_range_x() const;
  S half_range_y() const;

  /// Update the height map using the points, iterators and functors
  void resetHeights();
  void updateHeightsByPointCloud3D(
      const PointCloud& points,
      const Eigen::Isometry3f* points_to_heightmap_frame = nullptr);
  template <typename Point3DIterator>
  void updateHeightsByPoint3DIterator(
      const Point3DIterator& begin, const Point3DIterator& end,
      const Eigen::Isometry3f* points_to_heightmap_frame = nullptr);
  using PointGenerationFunc = typename FlatHeightMap<S>::PointGenerationFunc;
  void updateHeightsByPointGenerationFunctor(
      const PointGenerationFunc& point_generator, int n_points);

  /// Update the heights by visitor that updates the bottom layer
  using UpdateBottomLayerHeightsFunctor = std::function<bool(
      const Pixel& bottom_pixel, const Point2D<S>& box_bottom_center,
      uint16_t old_height_in_mm, uint16_t& new_height_in_mm)>;
  void updateHeightsByBottomLayerUpdateFunctor(
      const UpdateBottomLayerHeightsFunctor& visitor,
      const PixelSpaceROI* roi = nullptr);

 private:
  void rebuildNextLayer(const FlatHeightMap<S>& down_layer,
                        FlatHeightMap<S>& up_layer);
  void updateEveryLayersFromBottom();
};

}  // namespace heightmap
}  // namespace fcl

#include "fcl/geometry/heightmap/layered_heightmap-inl.h"