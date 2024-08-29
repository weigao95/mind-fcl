//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

#include <cstdint>

namespace fcl {
namespace octree2 {

/// Node index in a flatten vector, unsigned similar to size_t
using OctreeNodeIndex = std::uint32_t;
constexpr OctreeNodeIndex kInvalidNodeIndex =
    std::numeric_limits<OctreeNodeIndex>::max();

/// Indexing the child within a OctreeNode, in range [0, 7]
/// Thus, 8 is used as invalid number
using OctreeChildIndex = std::uint8_t;
constexpr OctreeChildIndex kNumberChildOfOctant = 8;
constexpr OctreeChildIndex kInvalidChildIndex = 8;

/// When a user-defined resolution is given (e.g., 2[mm]) to make a octree,
/// the world is discretized into a volume grid with that resolution. Then, a
/// parent level node consists of 8 children nodes which form a 2x2x2 cube.
///
/// To improve the efficient, the finest level is not directly stored in the
/// Octree. Actually, the LeafNode in the Octree represents a 2x2x2 grid, thus
/// its resolution is 2x user-provided. Then, other OctreeNode built upon this
/// LeafNode would be InnerNode.
///
/// This enum distinguish types of node at different resolutions.
enum class OctreeLayerNodeType : std::uint8_t {
  InnerNode = 0,
  LeafNode = 1,
  RealLeafNoNode = 2,
};

/// Octree is sparse representation of dense volumetric field. The voxel type
/// is actually from dense volume field. This is very similar to Pixel in
/// HeightMap, which can be regarded as an 1-channel depth image.
/// As Pixel type is exactly the pixel on that image in the usual image
/// coordinate:
///     x corresponds to rows and range from [0, n_rows), n_rows = map_shape_x
///     y corresponds to rows and range from [0, n_cols), n_rows = map_shape_y
/// The voxel is defined as a simple extension of the Pixel in HeightMap with
/// extension to include the z axis.
/// On the other hand, the Octree represents a 3D-coordinate geometry
/// consists of many boxes. z-axis is the height. The ranges is:
///    [-half_range_x, half_range_x] x [-h_y, h_y] x [-h_z, h_z]
/// where: half_range_x/y/z = bottom_half_shape * resolution_x/y/z.
/// Each voxel is a 3d Box in the world.
struct OctreeVoxel {
  std::array<std::uint16_t, 3> xyz{};
  std::uint16_t x() const { return xyz[0]; }
  std::uint16_t y() const { return xyz[1]; }
  std::uint16_t z() const { return xyz[2]; }
  std::uint16_t& x() { return xyz[0]; }
  std::uint16_t& y() { return xyz[1]; }
  std::uint16_t& z() { return xyz[2]; }
};

}  // namespace octree2
}  // namespace fcl
