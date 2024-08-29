#include <gtest/gtest.h>

#include "fcl/geometry/heightmap/layered_heightmap.h"
#include "fcl/geometry/shape/utility.h"

namespace fcl {
namespace heightmap {

template <typename S>
void generateRepresentativeHeightMapConfiguration(
    std::vector<LayeredHeightMap<S>>& height_map_to_test) {
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 256));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 512));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 0.002, 256, 256));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 0.002, 512, 512));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 0.001, 256, 512));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 0.001, 512, 256));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.002, 0.001, 256, 512));
  height_map_to_test.emplace_back(LayeredHeightMap<S>(0.001, 0.002, 512, 256));
}

template <typename S>
void updateByRandomPointCloud(std::size_t n_points, PointCloud& point_cloud,
                              LayeredHeightMap<S>& heightmap,
                              float& max_cloud_height,
                              bool exclude_out_range_points) {
  point_cloud.clear();
  point_cloud.reserve(n_points);
  max_cloud_height = -std::numeric_limits<float>::infinity();
  while (point_cloud.size() < n_points) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    // Remove the point if it is too far
    if (exclude_out_range_points &&
        (std::abs(point.x()) >= heightmap.bottom().half_range_x() ||
         std::abs(point.y()) >= heightmap.bottom().half_range_y())) {
      continue;
    }

    // Some of the points are below check xOy plane
    const float z_height = point[2] + 0.9;
    point_cloud.push_back(point[0], point[1], z_height);
    max_cloud_height = std::max(max_cloud_height, z_height);
  }

  // Perform update
  heightmap.updateHeightsByPointCloud3D(point_cloud);
}

template <typename S>
bool testAABBContainmentOfEachPixel(const LayeredHeightMap<S>& height_map,
                                    S expected_max_height_m) {
  // Compute the heightmap aabb
  AABB<S> global_aabb;
  global_aabb.min_ =
      Vector3<S>(-height_map.half_range_x(), -height_map.half_range_y(), 0);
  global_aabb.max_ =
      Vector3<S>(height_map.half_range_x(), height_map.half_range_y(),
                 expected_max_height_m);

  // For each box
  Pixel pixel_xy;
  AABB<S> pixel_aabb;
  Box<S> pixel_box;
  Transform3<S> pixel_box_tf;
  for (const auto& layer : height_map.layers()) {
    for (uint16_t y = 0; y < layer.full_shape_y(); y++) {
      for (uint16_t x = 0; x < layer.full_shape_x(); x++) {
        pixel_xy = Pixel{x, y};
        bool ok = layer.pixelToBox(pixel_xy, pixel_aabb);
        if (!ok) {
          EXPECT_TRUE(layer.pixelHeight(pixel_xy) == 0);
          EXPECT_FALSE(layer.pixelToBox(pixel_xy, pixel_box, pixel_box_tf));
        } else {
          if (!global_aabb.contain(pixel_aabb)) {
            return false;
          }

          // Compute the box and compare
          EXPECT_TRUE(layer.pixelToBox(pixel_xy, pixel_box, pixel_box_tf));
          {
            Box<S> box_from_aabb;
            Transform3<S> box_tf_from_aabb;
            constructBox(pixel_aabb, Transform3<S>::Identity(), box_from_aabb,
                         box_tf_from_aabb);
            if (std::is_same<S, double>::value) {
              EXPECT_NEAR((pixel_box.side - box_from_aabb.side).norm(), 0.0,
                          1e-7);
              EXPECT_NEAR(
                  (pixel_box_tf.translation() - box_tf_from_aabb.translation())
                      .norm(),
                  0.0, 1e-7);
            } else {
              EXPECT_NEAR((pixel_box.side - box_from_aabb.side).norm(), 0.0,
                          2e-7);
              EXPECT_NEAR(
                  (pixel_box_tf.translation() - box_tf_from_aabb.translation())
                      .norm(),
                  0.0, 1e-7);
            }
          }
        }
      }
    }
  }

  // Done
  return true;
}

template <typename S>
void testCoordinateConversion() {
  // Height map dims
  std::vector<LayeredHeightMap<S>> height_map_to_test;
  generateRepresentativeHeightMapConfiguration(height_map_to_test);

  // For test the height map basic info
  for (const auto& height_map : height_map_to_test) {
    EXPECT_TRUE(height_map.n_layers() >= 2);
    EXPECT_TRUE(height_map.top().half_shape_x() >= 1);
    EXPECT_TRUE(height_map.top().half_shape_y() >= 1);
    EXPECT_TRUE(height_map.top().half_shape_x() == 1 ||
                height_map.top().half_shape_y() == 1);
    EXPECT_NEAR(height_map.top().half_range_x(),
                height_map.bottom().half_range_x(), 1e-7);
    EXPECT_NEAR(height_map.top().half_range_y(),
                height_map.bottom().half_range_y(), 1e-7);
  }

  // Start testing by points
  constexpr int random_test_points = 1000000;
  for (const auto& height_map : height_map_to_test) {
    for (int test_i = 0; test_i < random_test_points; test_i++) {
      Point2D<S> point2d_init = Point2D<S>::Random();
      // Sometimes this would be exactly 1 from random
      Point2D<S> point2d = point2d_init * S(1.0 - 1e-6);
      point2d.x() *= height_map.half_range_x();
      point2d.y() *= height_map.half_range_y();

      // Ok for each layer
      Pixel point_pixel;
      for (const auto& layer : height_map.layers()) {
        bool ok = layer.point2DToPixel(point2d, point_pixel);
        EXPECT_TRUE(ok);
        EXPECT_TRUE(layer.isPixelInRange(point_pixel));
      }
    }
  }
}

template <typename S>
void testAABBContainment() {
  // Height map dims
  std::vector<LayeredHeightMap<S>> height_map_to_test;
  generateRepresentativeHeightMapConfiguration(height_map_to_test);

  // pixel.x + pixel.y heightmap
  uint16_t product_x_by = 1;
  auto update_heightmap_x_plus_y =
      [&product_x_by](const Pixel& pixel, const Point2D<S>& box_bottom_center,
                      uint16_t old_height_in_mm,
                      uint16_t& new_height_in_mm) -> bool {
    (void)(box_bottom_center);
    (void)(old_height_in_mm);
    new_height_in_mm = product_x_by * pixel.x + pixel.y;
    return false;
  };

  // Test loop
  for (auto& height_map : height_map_to_test) {
    height_map.resetHeights();
    product_x_by = 1;
    height_map.updateHeightsByBottomLayerUpdateFunctor(
        update_heightmap_x_plus_y);
    uint16_t expected_max_height_mm = height_map.bottom().full_shape_x() +
                                      height_map.bottom().full_shape_y() - 2;
    S expected_max_height_m = S(expected_max_height_mm) * S(0.001);
    EXPECT_EQ(expected_max_height_mm, height_map.height_upper_bound_mm());
    EXPECT_NEAR(expected_max_height_m, height_map.height_upper_bound_meter(),
                1e-6);
    EXPECT_TRUE(
        testAABBContainmentOfEachPixel<S>(height_map, expected_max_height_m));

    // Test again with taller value
    product_x_by = 2;
    height_map.updateHeightsByBottomLayerUpdateFunctor(
        update_heightmap_x_plus_y);
    uint16_t expected_max_height_mm_false =
        2 * height_map.bottom().full_shape_x() +
        height_map.bottom().full_shape_y() - 4;
    S expected_max_height_m_false = S(expected_max_height_mm_false) * S(0.001);
    EXPECT_FALSE(testAABBContainmentOfEachPixel<S>(
        height_map, expected_max_height_m_false));
    expected_max_height_mm = 2 * height_map.bottom().full_shape_x() +
                             height_map.bottom().full_shape_y() - 3;
    expected_max_height_m = S(expected_max_height_mm) * S(0.001);
    EXPECT_EQ(expected_max_height_mm, height_map.height_upper_bound_mm());
    EXPECT_NEAR(expected_max_height_m, height_map.height_upper_bound_meter(),
                1e-6);
    EXPECT_TRUE(
        testAABBContainmentOfEachPixel<S>(height_map, expected_max_height_m));
  }
}

}  // namespace heightmap
}  // namespace fcl

GTEST_TEST(LayeredHeightMapTest, BasicCoordinateConversionTest) {
  fcl::heightmap::testCoordinateConversion<float>();
  fcl::heightmap::testCoordinateConversion<double>();
}

GTEST_TEST(LayeredHeightMapTest, AABB_ContainmentTest) {
  fcl::heightmap::testAABBContainment<float>();
  fcl::heightmap::testAABBContainment<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}