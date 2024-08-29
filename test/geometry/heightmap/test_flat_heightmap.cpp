#include <gtest/gtest.h>
#include "fcl/geometry/shape/utility.h"
#include "fcl/geometry/heightmap/flat_heightmap.h"

namespace fcl {
namespace heightmap {

template <typename S>
void generateRepresentativeHeightMapConfiguration(
    std::vector<FlatHeightMap<S>>& height_map_to_test) {
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 256));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 512));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 0.002, 256, 256));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 0.002, 512, 512));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 0.001, 256, 512));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 0.001, 512, 256));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.002, 0.001, 256, 512));
  height_map_to_test.emplace_back(FlatHeightMap<S>(0.001, 0.002, 512, 256));
}

template <typename S>
void updateByRandomPointCloud(std::size_t n_points, PointCloud& point_cloud,
                              FlatHeightMap<S>& heightmap,
                              float& max_cloud_height,
                              bool exclude_out_range_points) {
  point_cloud.clear();
  point_cloud.reserve(n_points);
  max_cloud_height = -std::numeric_limits<float>::infinity();
  while (point_cloud.size() < n_points) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    // Remove the point if it is too far
    if (exclude_out_range_points &&
        (std::abs(point.x()) >= heightmap.half_range_x() ||
         std::abs(point.y()) >= heightmap.half_range_y())) {
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
bool testAABBContainmentOfEachPixel(const FlatHeightMap<S>& height_map,
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
  for (uint16_t y = 0; y < height_map.full_shape_y(); y++) {
    for (uint16_t x = 0; x < height_map.full_shape_x(); x++) {
      pixel_xy = Pixel{x, y};
      bool ok = height_map.pixelToBox(pixel_xy, pixel_aabb);
      if (!ok) {
        EXPECT_TRUE(height_map.pixelHeight(pixel_xy) == 0);
        EXPECT_FALSE(height_map.pixelToBox(pixel_xy, pixel_box, pixel_box_tf));
      } else {
        if (!global_aabb.contain(pixel_aabb)) {
          return false;
        }

        // Compute the box and compare
        EXPECT_TRUE(height_map.pixelToBox(pixel_xy, pixel_box, pixel_box_tf));
        {
          Box<S> box_from_aabb;
          Transform3<S> box_tf_from_aabb;
          constructBox(pixel_aabb, Transform3<S>::Identity(), box_from_aabb,
                       box_tf_from_aabb);
          if (std::is_same<S, double>::value) {
            EXPECT_NEAR((pixel_box.side - box_from_aabb.side).norm(), 0.0, 1e-7);
            EXPECT_NEAR((pixel_box_tf.translation() - box_tf_from_aabb.translation()).norm(), 0.0, 1e-7);
          } else {
            EXPECT_NEAR((pixel_box.side - box_from_aabb.side).norm(), 0.0, 2e-7);
            EXPECT_NEAR((pixel_box_tf.translation() - box_tf_from_aabb.translation()).norm(), 0.0, 1e-7);
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
  std::vector<FlatHeightMap<S>> height_map_to_test;
  generateRepresentativeHeightMapConfiguration(height_map_to_test);

  // Start testing
  constexpr int random_test_points = 1000000;
  for (const auto& height_map : height_map_to_test) {
    for (int test_i = 0; test_i < random_test_points; test_i++) {
      Point2D<S> point2d_init = Point2D<S>::Random();
      // Sometimes this would be exactly 1 from random
      Point2D<S> point2d = point2d_init * S(1.0 - 1e-6);
      point2d.x() *= height_map.half_range_x();
      point2d.y() *= height_map.half_range_y();
      Pixel point_pixel;
      bool ok = height_map.point2DToPixel(point2d, point_pixel);
      EXPECT_TRUE(ok);
      EXPECT_TRUE(height_map.isPixelInRange(point_pixel));

      // To point and check diff
      Point2D<S> point2d_center, point2d_topleft, point2d_bottomright;
      ok = height_map.pixelToPoint2D(point_pixel, point2d_center,
                                     PixelToPoint2DType::Center);
      EXPECT_TRUE(ok);
      ok = height_map.pixelToPoint2D(point_pixel, point2d_topleft,
                                     PixelToPoint2DType::TopLeft);
      EXPECT_TRUE(ok);
      ok = height_map.pixelToPoint2D(point_pixel, point2d_bottomright,
                                     PixelToPoint2DType::BottomRight);
      EXPECT_TRUE(ok);

      // Check the difference between point2d
      // clang-format off
      EXPECT_TRUE(std::abs(point2d.x() - point2d_center.x()) <= 0.5 * height_map.resolution_x() + 1e-6);
      EXPECT_TRUE(std::abs(point2d.y() - point2d_center.y()) <= 0.5 * height_map.resolution_y() + 1e-6);
      EXPECT_TRUE(std::abs(point2d.x() - point2d_topleft.x()) <= height_map.resolution_x() + 1e-6);
      EXPECT_TRUE(std::abs(point2d.y() - point2d_topleft.y()) <= height_map.resolution_y() + 1e-6);
      EXPECT_TRUE(std::abs(point2d.x() - point2d_bottomright.x()) <= height_map.resolution_x() + 1e-6);
      EXPECT_TRUE(std::abs(point2d.y() - point2d_bottomright.y()) <= height_map.resolution_y() + 1e-6);
      // clang-format on
    }
  }
}

template <typename S>
void testAABBContainment() {
  // Height map dims
  std::vector<FlatHeightMap<S>> height_map_to_test;
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
    height_map.updateHeightsByFunctor(update_heightmap_x_plus_y);
    uint16_t expected_max_height_mm =
        height_map.full_shape_x() + height_map.full_shape_y() - 2;
    S expected_max_height_m = S(expected_max_height_mm) * S(0.001);
    EXPECT_EQ(expected_max_height_mm, height_map.height_upper_bound_mm());
    EXPECT_NEAR(expected_max_height_m, height_map.height_upper_bound_meter(),
                1e-6);
    EXPECT_TRUE(
        testAABBContainmentOfEachPixel<S>(height_map, expected_max_height_m));

    // Test again with taller value
    product_x_by = 2;
    height_map.updateHeightsByFunctor(update_heightmap_x_plus_y);
    uint16_t expected_max_height_mm_false =
        2 * height_map.full_shape_x() + height_map.full_shape_y() - 4;
    S expected_max_height_m_false = S(expected_max_height_mm_false) * S(0.001);
    EXPECT_FALSE(testAABBContainmentOfEachPixel<S>(
        height_map, expected_max_height_m_false));
    expected_max_height_mm =
        2 * height_map.full_shape_x() + height_map.full_shape_y() - 3;
    expected_max_height_m = S(expected_max_height_mm) * S(0.001);
    EXPECT_EQ(expected_max_height_mm, height_map.height_upper_bound_mm());
    EXPECT_NEAR(expected_max_height_m, height_map.height_upper_bound_meter(),
                1e-6);
    EXPECT_TRUE(
        testAABBContainmentOfEachPixel<S>(height_map, expected_max_height_m));
  }
}

template <typename S>
void testRandomPointCloudAABBContainment() {
  // Height map dims
  std::vector<FlatHeightMap<S>> height_map_to_test;
  generateRepresentativeHeightMapConfiguration(height_map_to_test);

  // Start testing
  constexpr int random_test_points = 1000000;
  PointCloud point_cloud;
  float max_cloud_height;
  for (auto& height_map : height_map_to_test) {
    // All points within the range
    height_map.resetHeights();
    updateByRandomPointCloud(random_test_points, point_cloud, height_map,
                             max_cloud_height, true);
    EXPECT_TRUE(testAABBContainmentOfEachPixel<S>(
        height_map, height_map.height_upper_bound_meter()));
    S max_height_from_heightmap = height_map.height_upper_bound_meter();
    EXPECT_NEAR(max_cloud_height, max_height_from_heightmap, 0.001);

    // Not all points within the range
    height_map.resetHeights();
    updateByRandomPointCloud(random_test_points, point_cloud, height_map,
                             max_cloud_height, false);
    EXPECT_TRUE(testAABBContainmentOfEachPixel<S>(
        height_map, height_map.height_upper_bound_meter()));
    max_height_from_heightmap = height_map.height_upper_bound_meter();
    EXPECT_TRUE(max_cloud_height >= max_height_from_heightmap - 0.001);
  }
}

template <typename S>
void testRegionOfInterest() {
  // Height map dims
  std::vector<FlatHeightMap<S>> height_map_to_test;
  generateRepresentativeHeightMapConfiguration(height_map_to_test);

  auto update_heightmap_x_plus_y =
      [](const Pixel& pixel, const Point2D<S>& box_bottom_center,
         uint16_t old_height_in_mm, uint16_t& new_height_in_mm) -> bool {
    (void)(box_bottom_center);
    (void)(old_height_in_mm);
    new_height_in_mm = pixel.x + pixel.y;
    return false;
  };

  bool inspect_ok = true;
  auto inspect_heightmap_x_plus_y =
      [&inspect_ok](const Pixel& pixel, const Point2D<S>& box_bottom_center,
                    uint16_t height_in_mm) -> bool {
    (void)(box_bottom_center);
    if (height_in_mm != pixel.x + pixel.y) {
      inspect_ok = false;
      return true;
    } else {
      return false;
    }
  };

  // Start testing
  constexpr int random_test_points = 1000000;
  PointCloud point_cloud;
  float max_cloud_height;
  for (auto& height_map : height_map_to_test) {
    // All points within the range
    height_map.resetHeights();
    updateByRandomPointCloud(random_test_points, point_cloud, height_map,
                             max_cloud_height, true);

    // Overwrite out some of them
    PixelSpaceROI roi;
    roi.top_left.x = uint16_t(0.1 * height_map.full_shape_x());
    roi.top_left.y = uint16_t(0.1 * height_map.full_shape_y());
    roi.bottom_right.x = uint16_t(0.9 * height_map.full_shape_x());
    roi.bottom_right.y = uint16_t(0.9 * height_map.full_shape_y());
    height_map.updateHeightsByFunctor(update_heightmap_x_plus_y, &roi);

    // Check it
    inspect_ok = true;
    height_map.visitHeightMap(inspect_heightmap_x_plus_y, &roi);
    EXPECT_TRUE(inspect_ok);

    // Larger roi, should be false
    inspect_ok = true;
    height_map.visitHeightMap(inspect_heightmap_x_plus_y);
    EXPECT_FALSE(inspect_ok);

    // Try containment
    auto containment_ok = testAABBContainmentOfEachPixel(
        height_map, height_map.height_upper_bound_meter());
    EXPECT_TRUE(containment_ok);
  }
}

}  // namespace heightmap
}  // namespace fcl

GTEST_TEST(HeightMapTest, BasicCoordinateConversionTest) {
  fcl::heightmap::testCoordinateConversion<float>();
  fcl::heightmap::testCoordinateConversion<double>();
}

GTEST_TEST(HeightMapTest, AABB_ContainmentTest) {
  fcl::heightmap::testAABBContainment<float>();
  fcl::heightmap::testAABBContainment<double>();
}

GTEST_TEST(HeightMapTest, RandomCloudAABB_Containment) {
  fcl::heightmap::testRandomPointCloudAABBContainment<float>();
  fcl::heightmap::testRandomPointCloudAABBContainment<double>();
}

GTEST_TEST(HeightMapTest, ROI_Test_0) {
  fcl::heightmap::testRegionOfInterest<float>();
  fcl::heightmap::testRegionOfInterest<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}