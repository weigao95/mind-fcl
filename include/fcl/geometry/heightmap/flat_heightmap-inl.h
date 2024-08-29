#pragma once

namespace fcl {
namespace heightmap {

/// Constructor and build from points
template <typename S>
FlatHeightMap<S>::FlatHeightMap(S resolution_x, S resolution_y,
                                uint16_t half_map_shape_x,
                                uint16_t half_map_shape_y)
    : xy_info{{resolution_x, half_map_shape_x, uint16_t(half_map_shape_x * 2)},
              {resolution_y, half_map_shape_y, uint16_t(half_map_shape_y * 2)}},
      height_upper_bound_in_mm(0) {
  const std::size_t map_flat_size =
      xy_info[0].full_shape * xy_info[1].full_shape;
  heightmap_in_mm.resize(map_flat_size, 0);
}

template <typename S>
FlatHeightMap<S>::FlatHeightMap(S resolution, uint16_t half_map_shape)
    : FlatHeightMap(resolution, resolution, half_map_shape, half_map_shape) {}

/// Pixel and point2D query
template <typename S>
uint16_t FlatHeightMap<S>::pixelHeight(const Pixel& pixel) const {
  // Do NOT check zero for unsigned
  if (!isPixelInRange(pixel)) {
    return 0;
  }

  // OK to obtain height
  const std::size_t flat_index = pixelToFlatIndex(pixel);
  return heightmap_in_mm[flat_index];
}

template <typename S>
bool FlatHeightMap<S>::isPixelInRange(const Pixel& pixel) const {
  // Check range, do not check zero for unsigned type
  return (pixel.x < full_shape_x() && pixel.y < full_shape_y());
}

template <typename S>
bool FlatHeightMap<S>::point2DToPixel(const Point2D<S>& point,
                                      Pixel& pixel) const {
  const int x = floor(point[0] / resolution_x()) + half_shape_x();
  const int y = floor(point[1] / resolution_y()) + half_shape_y();
  bool in_range =
      (x >= 0 && x < full_shape_x() && y >= 0 && y < full_shape_y());
  pixel.x = static_cast<uint16_t>(x);
  pixel.y = static_cast<uint16_t>(y);
  return in_range;
}

template <typename S>
bool FlatHeightMap<S>::pixelToPoint2D(const Pixel& pixel, Point2D<S>& point,
                                      PixelToPoint2DType query_type) const {
  // Compute the point despite potential out-of-range
  pixelToPoint2DUnchecked(pixel, point, query_type);
  return isPixelInRange(pixel);
}

template <typename S>
bool FlatHeightMap<S>::aabbToPixelROI(const AABB<S>& aabb,
                                      PixelSpaceROI& roi) const {
  auto clamp = [](int x, int ub) -> int {
    if (x <= 0)
      return 0;
    else if (x >= ub)
      return ub;
    return x;
  };
  auto project_point2D = [this, &clamp](const Point2D<S>& point,
                                        Pixel& pixel) -> bool {
    const int x = floor(point[0] / resolution_x()) + half_shape_x();
    const int y = floor(point[1] / resolution_y()) + half_shape_y();
    pixel.x = static_cast<uint16_t>(clamp(x, full_shape_x() - 1));
    pixel.y = static_cast<uint16_t>(clamp(y, full_shape_y() - 1));
    bool is_projected = pixel.x != x || pixel.y != y;
    return is_projected;
  };

  assert(aabb.min_.x() <= aabb.max_.x());
  assert(aabb.min_.y() <= aabb.max_.y());
  bool top_left_is_projected =
      project_point2D({aabb.min_.x(), aabb.min_.y()}, roi.top_left);
  bool bottom_right_is_projected =
      project_point2D({aabb.max_.x(), aabb.max_.y()}, roi.bottom_right);
  if (top_left_is_projected && bottom_right_is_projected) {
    bool no_overlap = (roi.top_left.x == roi.bottom_right.x) ||
                      (roi.top_left.y == roi.bottom_right.y);
    return !no_overlap;
  }
  return true;
}

template <typename S>
void FlatHeightMap<S>::pixelToPoint2DUnchecked(
    const Pixel& pixel, Point2D<S>& point,
    PixelToPoint2DType query_type) const {
  // Compute the point despite potential out-of-range
  switch (query_type) {
    case PixelToPoint2DType::TopLeft: {
      point = {(pixel.x - half_shape_x()) * resolution_x(),
               (pixel.y - half_shape_y()) * resolution_y()};
      break;
    }
    case PixelToPoint2DType::BottomRight: {
      point = {(pixel.x - half_shape_x() + 1.0) * resolution_x(),
               (pixel.y - half_shape_y() + 1.0) * resolution_y()};
      break;
    }
    case PixelToPoint2DType::Center:
    default:
      point = {(pixel.x - half_shape_x() + 0.5) * resolution_x(),
               (pixel.y - half_shape_y() + 0.5) * resolution_y()};
  }
}

template <typename S>
std::size_t FlatHeightMap<S>::pixelToFlatIndex(const Pixel& pixel) const {
  return pixel.y * full_shape_x() + pixel.x;
}

template <typename S>
bool FlatHeightMap<S>::pixelToBox(const Pixel& pixel, Box<S>& box,
                                  Vector3<S>& box_center) const {
  // Check range
  if (!isPixelInRange(pixel)) {
    return false;
  }

  // The box is defined by FULL size (not half)
  // But we won't return a box if height is zero
  const uint16_t height_in_mm = pixelHeight(pixel);
  if (height_in_mm == 0) {
    return false;
  }

  const S height_in_meter = static_cast<S>(height_in_mm) * S(0.001);
  box = Box<S>(resolution_x(), resolution_y(), height_in_meter);

  // Compute the offset
  Point2D<S> bottom_center;
  pixelToPoint2DUnchecked(pixel, bottom_center, PixelToPoint2DType::Center);
  box_center.x() = bottom_center.x();
  box_center.y() = bottom_center.y();
  box_center.z() = S(0.5) * height_in_meter;
  return true;
}

template <typename S>
bool FlatHeightMap<S>::pixelToBox(
    const Pixel& pixel, Box<S>& box,
    Transform3<S>& box_tf_in_heightmap_frame) const {
  box_tf_in_heightmap_frame.linear().setIdentity();
  Vector3<S> box_center;
  auto ok = pixelToBox(pixel, box, box_center);
  box_tf_in_heightmap_frame.translation() = box_center;
  return ok;
}

template <typename S>
bool FlatHeightMap<S>::pixelToBox(const Pixel& pixel, AABB<S>& aabb) const {
  // Check range
  if (!isPixelInRange(pixel)) {
    return false;
  }

  // Height
  const uint16_t height_in_mm = pixelHeight(pixel);
  if (height_in_mm == 0) {
    return false;
  }

  // Valid box
  const S height_in_meter = static_cast<S>(height_in_mm) * S(0.001);
  aabb.min_.z() = S(0.0);
  aabb.max_.z() = height_in_meter;

  // x/y
  Point2D<S> bottom_center;
  pixelToPoint2DUnchecked(pixel, bottom_center, PixelToPoint2DType::Center);
  aabb.min_.x() = bottom_center.x() - S(0.5) * resolution_x();
  aabb.max_.x() = bottom_center.x() + S(0.5) * resolution_x();
  aabb.min_.y() = bottom_center.y() - S(0.5) * resolution_y();
  aabb.max_.y() = bottom_center.y() + S(0.5) * resolution_y();
  return true;
}

/// Update the height map using the points, iterators or functors
template <typename S>
void FlatHeightMap<S>::resetHeights() {
  std::fill(heightmap_in_mm.begin(), heightmap_in_mm.end(), 0);
  height_upper_bound_in_mm = 0;
}

template <typename S>
void FlatHeightMap<S>::updateHeightsByPointCloud3D(
    const PointCloud& points,
    const Eigen::Isometry3f* points_to_heightmap_frame) {
  auto generation_fn = [&](int index, S& x, S& y, S& z) -> void {
    if (points_to_heightmap_frame == nullptr) {
      x = points[index].x();
      y = points[index].y();
      z = points[index].z();
    } else {
      Eigen::Vector3f point3f{points[index].x(), points[index].x(),
                              points[index].x()};
      Eigen::Vector3f point_tf = (*points_to_heightmap_frame) * (point3f);
      x = point_tf.x();
      y = point_tf.y();
      z = point_tf.z();
    }
  };
  updateHeightsByPointGenerationFunctor(generation_fn, points.size());
  /*updateHeightsByPoint3DIterator(points.begin(), points.end(),
                                 points_to_heightmap_frame);*/
}

template <typename S>
template <typename Point3DIterator>
void FlatHeightMap<S>::updateHeightsByPoint3DIterator(
    const Point3DIterator& begin, const Point3DIterator& end,
    const Eigen::Isometry3f* points_to_heightmap_frame) {
  for (auto iter = begin; iter != end; ++iter) {
    // Compute point
    Eigen::Vector3f point3d{iter->x(), iter->y(), iter->z()};
    if (points_to_heightmap_frame != nullptr) {
      point3d = (*points_to_heightmap_frame) * point3d;
    }

    // Do not handle negative z
    if (point3d.z() < 0) continue;

    // Map to image
    Pixel pixel;
    bool in_range = point2DToPixel(Point2D<S>(point3d.x(), point3d.y()), pixel);
    if (in_range) {
      const auto index = pixel.y * full_shape_x() + pixel.x;
      const int z = static_cast<uint16_t>(point3d.z() * 1000);
      if (heightmap_in_mm[index] < z)
        heightmap_in_mm[index] = static_cast<uint16_t>(z);
      if (height_upper_bound_in_mm < heightmap_in_mm[index])
        height_upper_bound_in_mm = heightmap_in_mm[index];
    }
  }
}

template <typename S>
void FlatHeightMap<S>::updateHeightsByPointGenerationFunctor(
    const PointGenerationFunc& point_generator, int n_points) {
  for (auto i = 0; i < n_points; i++) {
    S point_x, point_y, point_z;
    point_generator(i, point_x, point_y, point_z);

    // Do not handle negative z
    if (point_z < 0) continue;

    // Map to image
    Pixel pixel;
    bool in_range = point2DToPixel(Point2D<S>(point_x, point_y), pixel);
    if (in_range) {
      const auto index = pixel.y * full_shape_x() + pixel.x;
      const int z = static_cast<uint16_t>(point_z * 1000);
      if (heightmap_in_mm[index] < z)
        heightmap_in_mm[index] = static_cast<uint16_t>(z);
      if (height_upper_bound_in_mm < heightmap_in_mm[index])
        height_upper_bound_in_mm = heightmap_in_mm[index];
    }
  }
}

template <typename S>
void FlatHeightMap<S>::visitHeightMap(const VisitHeightMapFunctor& visitor,
                                      const PixelSpaceROI* roi) const {
  // Obtain the range
  Pixel top_left, bottom_right;
  rectifyHeightMapROI(top_left, bottom_right, roi);

  // The visit loop without mutation
  for (uint16_t y = top_left.y; y <= bottom_right.y; y++) {
    for (uint16_t x = top_left.x; x <= bottom_right.x; x++) {
      Pixel pixel(x, y);
      Point2D<S> box_bottom_center;
      pixelToPoint2DUnchecked(pixel, box_bottom_center,
                              PixelToPoint2DType::Center);
      // Do not need checking here
      uint16_t height_in_mm = heightmap_in_mm[pixelToFlatIndex(pixel)];
      bool done = visitor(pixel, box_bottom_center, height_in_mm);
      if (done) return;
    }
  }
}

template <typename S>
PointCloud FlatHeightMap<S>::toPointCloud() const {
  octomap::Pointcloud cloud;
  cloud.reserve(heightmap_in_mm.size());
  auto append_point = [&cloud](const Pixel& pixel,
                               const Point2D<S>& box_bottom_center,
                               uint16_t height_in_mm) -> bool {
    (void)(pixel);
    cloud.push_back(point3d{float(box_bottom_center.x()),
                            float(box_bottom_center.y()),
                            float(height_in_mm) * float(0.001)});
    return false;
  };
  visitHeightMap(append_point);
  return cloud;
}

template <typename S>
void FlatHeightMap<S>::updateHeightsByFunctor(
    const UpdateHeightMapFunctor& visitor, const PixelSpaceROI* roi) {
  // Obtain the range
  Pixel top_left, bottom_right;
  rectifyHeightMapROI(top_left, bottom_right, roi);

  // The visit loop WITH mutation
  Pixel pixel;
  Point2D<S> box_bottom_center;
  uint16_t old_height_in_mm = 0;
  uint16_t new_height_in_mm = 0;
  bool done = false;

  for (uint16_t y = top_left.y; y <= bottom_right.y; y++) {
    for (uint16_t x = top_left.x; x <= bottom_right.x; x++) {
      pixel.x = x;
      pixel.y = y;
      pixelToPoint2DUnchecked(pixel, box_bottom_center,
                              PixelToPoint2DType::Center);
      // Do not need checking here
      old_height_in_mm = heightmap_in_mm[pixelToFlatIndex(pixel)];
      new_height_in_mm = old_height_in_mm;
      done =
          visitor(pixel, box_bottom_center, old_height_in_mm, new_height_in_mm);

      // Do need update the map
      if (old_height_in_mm != new_height_in_mm) {
        heightmap_in_mm[pixelToFlatIndex(pixel)] = new_height_in_mm;
        if (new_height_in_mm > height_upper_bound_in_mm) {
          height_upper_bound_in_mm = new_height_in_mm;
        }
      }

      // No-further visit
      if (done) {
        return;
      }
    }
  }
}

template <typename S>
S FlatHeightMap<S>::height_upper_bound_meter() const {
  return S(height_upper_bound_in_mm) * S(0.001);
}

template <typename S>
S FlatHeightMap<S>::half_range_x() const {
  return xy_info[0].resolution * xy_info[0].half_shape;
}

template <typename S>
S FlatHeightMap<S>::half_range_y() const {
  return xy_info[1].resolution * xy_info[1].half_shape;
}

template <typename S>
void FlatHeightMap<S>::rectifyHeightMapROI(Pixel& top_left, Pixel& bottom_right,
                                           const PixelSpaceROI* roi) const {
  // Note the containment semantic
  top_left = Pixel(0, 0);
  bottom_right = Pixel(full_shape_x() - 1, full_shape_y() - 1);
  if (roi != nullptr) {
    // Update with new roi
    top_left.x = std::max(roi->top_left.x, top_left.x);
    top_left.y = std::max(roi->top_left.y, top_left.y);
    bottom_right.x = std::min(roi->bottom_right.x, bottom_right.x);
    bottom_right.y = std::min(roi->bottom_right.y, bottom_right.y);

    // top_left <= bottom_right
    top_left.x = std::min(top_left.x, bottom_right.x);
    top_left.y = std::min(top_left.y, bottom_right.y);
  }
}

}  // namespace heightmap
}  // namespace fcl
