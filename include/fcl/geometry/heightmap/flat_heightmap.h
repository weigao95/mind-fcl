#pragma once
#include <array>
#include <memory>

#include "fcl/geometry/heightmap/heightmap_types.h"

namespace fcl {
namespace heightmap {

template <typename S>
class LayeredHeightMap;

/// HeightMap represents a 3D geometry consists of many boxes. These
/// boxes share a common bottom plane of the xOy, and the (positive) z-axis is
/// the height. The ranges in the xOy plane is:
///    [-half_range_x, half_range_x] x [-half_range_y, half_range_y]
/// where: half_range_x/y = map_half_shape_x/y * resolution_x/y.
/// Taking a top-down view from the z-axis, the Heightmap can be regarded as
/// an 1-channel image of depth. The 3D x-axis corresponds to row, which is also
/// the x-axis in usual image coordinate. Similar holds for y-axis and cols.
///
/// "Flat" heightmap implies this heightmap is stored as an flat "image".
/// Alternatively, we might consider a sparse or hierarchical representation.
/// \tparam S
template <typename S>
class FlatHeightMap {
 private:
  struct {
    S resolution;
    uint16_t half_shape;
    uint16_t full_shape;
  } xy_info[2];

  /// Vector with maximum height values for every bin (height map)
  /// Column-major representation (cols are continuous)
  std::vector<uint16_t> heightmap_in_mm;

  /// A upper bound of the heights in the map
  /// When there is no remove (or reduction of heights in update), this
  /// is also the max height. However, it becomes expensive to maintain
  /// max height during removal, and a upper_bound should be good enough
  /// for a global AABB of the heightmap.
  uint16_t height_upper_bound_in_mm{0};

  // Access of raw data
  friend class LayeredHeightMap<S>;

 public:
  FlatHeightMap(S resolution_x, S resolution_y, uint16_t half_map_shape_x,
                uint16_t half_map_shape_y);
  FlatHeightMap(S resolution, uint16_t half_map_shape);
  // Default destruction, copy/assign/move

  /// Pixel and point2D query
  // Obtain the height in [mm] of a pixel, return 0 if pixel out-of-range.
  uint16_t pixelHeight(const Pixel& pixel) const;
  bool isPixelInRange(const Pixel& pixel) const;

  // Compute the box of a given pixel
  // False would be return if height at pixel is zero
  bool pixelToBox(const Pixel& pixel, Box<S>& box,
                  Vector3<S>& box_center) const;
  bool pixelToBox(const Pixel& pixel, Box<S>& box,
                  Transform3<S>& box_tf_in_heightmap_frame) const;
  bool pixelToBox(const Pixel& pixel, AABB<S>& aabb) const;

  // Convert between point2D and pixel. Return value implies whether the input
  // is in the required range (for either point2D of pixel).
  bool point2DToPixel(const Point2D<S>& point, Pixel& pixel) const;
  bool pixelToPoint2D(
      const Pixel& pixel, Point2D<S>& point_output,
      PixelToPoint2DType query_type = PixelToPoint2DType::Center) const;

  /// Compute the region_of_interest of the aabb projected onto the heightmap.
  /// Return false indicates the area of region_of_interest is 0.
  bool aabbToPixelROI(const AABB<S>& aabb, PixelSpaceROI& roi) const;

 private:
  void pixelToPoint2DUnchecked(
      const Pixel& pixel, Point2D<S>& point_output,
      PixelToPoint2DType query_type = PixelToPoint2DType::Center) const;
  std::size_t pixelToFlatIndex(const Pixel& pixel) const;

  /// Update the height map using the points, iterators and functor
 public:
  void resetHeights();
  void updateHeightsByPointCloud3D(
      const PointCloud& points,
      const Eigen::Isometry3f* points_to_heightmap_frame = nullptr);
  template <typename Point3DIterator>
  void updateHeightsByPoint3DIterator(
      const Point3DIterator& begin, const Point3DIterator& end,
      const Eigen::Isometry3f* points_to_heightmap_frame = nullptr);
  using PointGenerationFunc = std::function<void(int index, S& x, S& y, S& z)>;
  void updateHeightsByPointGenerationFunctor(
      const PointGenerationFunc& point_generator, int n_points);

  /// Visit or update the height map using arbitrary functor
  /// The visitor/updater functor return a bool indicate
  /// WHETHER the visit CAN STOP
  /// visitor return TRUE implies NO further visit is necessary and we can
  /// stop now, return FALSE implies visit should continue.
  using VisitHeightMapFunctor = std::function<bool(
      const Pixel& pixel, const Point2D<S>& box_bottom_center,
      uint16_t height_in_mm)>;
  using UpdateHeightMapFunctor = std::function<bool(
      const Pixel& pixel, const Point2D<S>& box_bottom_center,
      uint16_t old_height_in_mm, uint16_t& new_height_in_mm)>;

  // Visit or update the height map
  // The visit is potentially limited by a region_of_interest
  void visitHeightMap(const VisitHeightMapFunctor& visitor,
                      const PixelSpaceROI* roi = nullptr) const;
  void updateHeightsByFunctor(const UpdateHeightMapFunctor& visitor,
                              const PixelSpaceROI* roi = nullptr);
  PointCloud toPointCloud() const;

  /// Simple access
 public:
  uint16_t half_shape_x() const { return xy_info[0].half_shape; }
  uint16_t half_shape_y() const { return xy_info[1].half_shape; }
  uint16_t full_shape_x() const { return xy_info[0].full_shape; }
  uint16_t full_shape_y() const { return xy_info[1].full_shape; }
  S resolution_x() const { return xy_info[0].resolution; }
  S resolution_y() const { return xy_info[1].resolution; }
  uint16_t height_upper_bound_mm() const { return height_upper_bound_in_mm; }
  S height_upper_bound_meter() const;
  S half_range_x() const;
  S half_range_y() const;

 private:
  void rectifyHeightMapROI(Pixel& top_left, Pixel& bottom_right,
                           const PixelSpaceROI* roi = nullptr) const;
};

}  // namespace heightmap
}  // namespace fcl

#include "fcl/geometry/heightmap/flat_heightmap-inl.h"