#pragma once

#include "fcl/geometry/octree/octomap/Pointcloud.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/math/bv/AABB.h"

namespace fcl {
namespace heightmap {

/// Borrow the point cloud type from octomap
using PointCloud = ::fcl::octomap::Pointcloud;
using point3d = ::fcl::octomap::point3d;

/// Heightmap can be regarded as an 1-channel depth image. Pixel type is
/// exactly the pixel on that image in the usual image coordinate:
///     x corresponds to rows and range from [0, n_rows), n_rows = map_shape_x
///     y corresponds to rows and range from [0, n_cols), n_rows = map_shape_y
/// On the other hand, the heightmap represents a 3D-coordinate geometry
/// consists of many boxes. z-axis is the height. The ranges in xOy plane is:
///    [-half_range_x, half_range_x] x [-half_range_y, half_range_y]
/// where: half_range_x/y = map_half_shape_x/y * resolution_x/y.
///
/// Each pixel is a 2d Box in the xOy plane, and we are usually interested in
/// its 4 corner points and center point.
struct Pixel {
  uint16_t x;  // [0, map_shape_x)
  uint16_t y;  // [0, map_shape_y)
  uint16_t row() const { return x; }
  uint16_t col() const { return y; }
  explicit Pixel() : x(0), y(0) {}
  Pixel(uint16_t x_in, uint16_t y_in) : x(x_in), y(y_in) {}
};

/// As mentioned above, each pixel is a 2d Box in the xOy plane, and we
/// are usually interested in its 4 corner points and center point. This enum
/// flags which point we are interested in.
enum class PixelToPoint2DType { Center, TopLeft, BottomRight };

/// 2D point type, usually in xOy plane
template <typename S>
using Point2D = fcl::Vector2<S>;

/// RegionOfInterest type in pixel coordinate, using both-include convention:
///              [top, bottom] X [left, right]
/// Note that bottom and right are INCLUDED
struct PixelSpaceROI {
  Pixel top_left;
  Pixel bottom_right;
};

/// The upper 16 bits of int store Pixel.x, and the lower 16 bits of int store
/// Pixel.y.
inline int encodePixel(const Pixel& pixel) {
  int code = pixel.x;
  code <<= 16;
  code |= pixel.y;
  return code;
}

inline Pixel decodePixel(int code) {
  Pixel pixel;
  pixel.x = (code >> 16) & 0xFFFF;
  pixel.y = code & 0xFFFF;
  return pixel;
}

}  // namespace heightmap
}  // namespace fcl