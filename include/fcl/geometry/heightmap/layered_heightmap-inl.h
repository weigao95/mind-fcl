#pragma once

namespace fcl {
namespace heightmap {

template <typename S>
LayeredHeightMap<S>::LayeredHeightMap(S bottom_resolution_x,
                                      S bottom_resolution_y,
                                      uint16_t bottom_half_map_shape_x,
                                      uint16_t bottom_half_map_shape_y) {
  // TODO(wei): consider shall we throw or rectify it here?
  uint16_t n_x = bottom_half_map_shape_x;
  uint16_t n_y = bottom_half_map_shape_y;
  assert(n_x > 0 && ((n_x & (n_x - 1)) == 0));
  assert(n_y > 0 && ((n_y & (n_y - 1)) == 0));

  // As layers_[0] is the top, we must construct in reversed order
  // To make it reversed, compute the top level resolution and shape
  S resolution_x = bottom_resolution_x;
  S resolution_y = bottom_resolution_y;
  while (true) {
    if (n_x > 1 && n_y > 1) {
      // Update the size and resolution
      resolution_x *= 2;
      resolution_y *= 2;
      n_x = n_x / 2;
      n_y = n_y / 2;
    } else {
      break;
    }
  }
  assert(n_x == 1 || n_y == 1);

  // Build from the top to bottom
  while (true) {
    // Make a new layer
    layers_.emplace_back(
        FlatHeightMap<S>(resolution_x, resolution_y, n_x, n_y));
    if (n_x != bottom_half_map_shape_x) {
      resolution_x *= S(0.5);
      resolution_y *= S(0.5);
      n_x *= 2;
      n_y *= 2;
    } else {
      assert(n_y == bottom_half_map_shape_y);
      break;
    }
  }
}

template <typename S>
LayeredHeightMap<S>::LayeredHeightMap(S bottom_resolution,
                                      uint16_t bottom_half_map_shape)
    : LayeredHeightMap(bottom_resolution, bottom_resolution,
                       bottom_half_map_shape, bottom_half_map_shape) {}

template <typename S>
bool LayeredHeightMap<S>::is_bottom_layer(std::size_t layer_idx) const {
  return layer_idx == (layers_.size() - 1);
}

template <typename S>
S LayeredHeightMap<S>::height_upper_bound_meter() const {
  return layers_[kTopLayerIndex].height_upper_bound_meter();
}

template <typename S>
uint16_t LayeredHeightMap<S>::height_upper_bound_mm() const {
  return layers_[kTopLayerIndex].height_upper_bound_mm();
}

template <typename S>
S LayeredHeightMap<S>::half_range_x() const {
  return layers_[kTopLayerIndex].half_range_x();
}

template <typename S>
S LayeredHeightMap<S>::half_range_y() const {
  return layers_[kTopLayerIndex].half_range_y();
}

template <typename S>
void LayeredHeightMap<S>::rebuildNextLayer(const FlatHeightMap<S>& down_layer,
                                           FlatHeightMap<S>& up_layer) {
  assert(up_layer.half_shape_x() * 2 == down_layer.half_shape_x());
  assert(up_layer.half_shape_y() * 2 == down_layer.half_shape_y());
  up_layer.resetHeights();
  auto update_up_layer_height = [&up_layer](const Pixel& pixel_down,
                                            const Point2D<S>& box_bottom_center,
                                            uint16_t height_in_mm) -> bool {
    (void)(box_bottom_center);
    // A pixel in up layer corresponds an 2x2 grid in down layer
    Pixel pixel_up;
    pixel_up.x = pixel_down.x / 2;
    pixel_up.y = pixel_down.y / 2;
    auto flat_index_in_up = up_layer.pixelToFlatIndex(pixel_up);
    if (height_in_mm > up_layer.heightmap_in_mm[flat_index_in_up]) {
      up_layer.heightmap_in_mm[flat_index_in_up] = height_in_mm;
    }

    // Visit the entire down map
    return false;
  };
  down_layer.visitHeightMap(update_up_layer_height);
  up_layer.height_upper_bound_in_mm = down_layer.height_upper_bound_in_mm;
}

template <typename S>
void LayeredHeightMap<S>::updateEveryLayersFromBottom() {
  assert(layers_.size() >= 1);
  for (std::size_t i = 0; i < layers_.size() - 1; i++) {
    std::size_t down_i = layers_.size() - 1 - i;
    std::size_t up_i = down_i - 1;
    rebuildNextLayer(layers_[down_i], layers_[up_i]);
  }
}

template <typename S>
void LayeredHeightMap<S>::resetHeights() {
  for (auto& layer : layers_) {
    layer.resetHeights();
  }
}

template <typename S>
void LayeredHeightMap<S>::updateHeightsByPointCloud3D(
    const PointCloud& points,
    const Eigen::Isometry3f* points_to_heightmap_frame) {
  updateHeightsByPoint3DIterator(points.begin(), points.end(),
                                 points_to_heightmap_frame);
}

template <typename S>
template <typename Point3DIterator>
void LayeredHeightMap<S>::updateHeightsByPoint3DIterator(
    const Point3DIterator& begin, const Point3DIterator& end,
    const Eigen::Isometry3f* points_to_heightmap_frame) {
  bottom_mutable().updateHeightsByPoint3DIterator(begin, end,
                                                  points_to_heightmap_frame);
  updateEveryLayersFromBottom();
}

template <typename S>
void LayeredHeightMap<S>::updateHeightsByPointGenerationFunctor(
    const PointGenerationFunc& point_generator, int n_points) {
  bottom_mutable().updateHeightsByPointGenerationFunctor(point_generator,
                                                         n_points);
  updateEveryLayersFromBottom();
}

template <typename S>
void LayeredHeightMap<S>::updateHeightsByBottomLayerUpdateFunctor(
    const UpdateBottomLayerHeightsFunctor& visitor, const PixelSpaceROI* roi) {
  bottom_mutable().updateHeightsByFunctor(visitor, roi);
  updateEveryLayersFromBottom();
}

}  // namespace heightmap
}  // namespace fcl