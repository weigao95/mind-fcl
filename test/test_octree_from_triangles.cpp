//
// Created by mech-mind_gw on 12/25/2023.
//
#include <gtest/gtest.h>

#include "octree_from_triangles.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void triangleToVoxelInstanceTest() {
  constexpr int test_n = 100;
  constexpr S resolution = 0.01;

  // Test loop
  Vector3<S> a, b, c;
  a.setRandom();
  b.setRandom();
  c.setRandom();
  TriangleP<S> triangle(a, b, c);
  std::unordered_set<std::int64_t> encoded_voxel_set;
  for (auto i = 0; i < test_n; i++) {
    triangle.a.setRandom();
    triangle.b.setRandom();
    triangle.c.setRandom();
    encoded_voxel_set.clear();
    makeTriangleSoupSurfaceVoxel({triangle}, resolution, encoded_voxel_set);
    EXPECT_FALSE(encoded_voxel_set.empty());
  }
}

template <typename S>
void rectangleAsSegmentsTest() {
  /* x-axis is downward in screen, y-axis points to right
   * top-left is (0, 0)
   * -----
   * |   |
   * |  ||
   * | | |
   * ||  |
   * |   |
   * -----
   * bottom left is (1, 0)
   * */
  const S rectangle_edge_length = 1.0;
  const int n_segments = 1000;
  const S segment_length = rectangle_edge_length / n_segments;
  constexpr S resolution = 0.001;  // 1 [mm]

  // Randomly apply a transformation
  fcl::Transform3<S> rectangle_tf;
  std::array<S, 6> xyz_extent{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  fcl::test::generateRandomTransform(xyz_extent, rectangle_tf);

  std::vector<TriangleP<S>> triangles;
  for (auto i = 0; i < n_segments; i++) {
    const S y_i = i * segment_length;
    const S y_i_plus_1 = y_i + segment_length;
    {
      Vector3<S> a_top_left(0, y_i, 0);
      Vector3<S> b_bottom_left(1.0, y_i, 0);
      Vector3<S> c_top_right(0, y_i_plus_1, 0);
      TriangleP<S> triangle(rectangle_tf * a_top_left,
                            rectangle_tf * b_bottom_left,
                            rectangle_tf * c_top_right);
      triangles.emplace_back(std::move(triangle));
    }

    {
      Vector3<S> a_top_right(0, y_i_plus_1, 0);
      Vector3<S> b_bottom_left(1.0, y_i, 0);
      Vector3<S> c_bottom_right(1.00, y_i_plus_1, 0);
      TriangleP<S> triangle(rectangle_tf * a_top_right,
                            rectangle_tf * b_bottom_left,
                            rectangle_tf * c_bottom_right);
      triangles.emplace_back(std::move(triangle));
    }
  }

  std::unordered_set<std::int64_t> encoded_voxel_set;
  makeTriangleSoupSurfaceVoxel(triangles, resolution, encoded_voxel_set);
  EXPECT_FALSE(encoded_voxel_set.empty());
}

}  // namespace fcl

GTEST_TEST(OctreeFromTrianglesTest, SurfaceInstanceTest) {
  fcl::triangleToVoxelInstanceTest<float>();
  fcl::triangleToVoxelInstanceTest<double>();
}

GTEST_TEST(OctreeFromTrianglesTest, RectangleAsSegmentsTest) {
  fcl::rectangleAsSegmentsTest<float>();
  fcl::rectangleAsSegmentsTest<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}