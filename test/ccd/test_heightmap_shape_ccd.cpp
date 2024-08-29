//
// Created by Wei Gao on 2024/6/18.
//
#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/detail/ccd/heightmap_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

std::shared_ptr<heightmap::PointCloud> generateRandomPointCloud(
    std::size_t numPoints) {
  auto cloud = std::make_shared<fcl::octomap::Pointcloud>();
  cloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    cloud->push_back(point[0], point[1], point[2] + 1);
  }
  return cloud;
}

template <typename S, typename Shape>
void heightmapShapeTranslationalCollisionRandomTest(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const TranslationalDisplacement<S>& hm1_displacement, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  // Construct the solver
  fcl::detail::TranslationalDisplacementHeightMapSolver<S> solver;

  // Check the collision between the shape and the small boxes that
  // make up the heightmap
  {
    ContinuousCollisionRequest<S> request;
    request.num_max_contacts = 100000;
    ContinuousCollisionResult<S> result;
    solver.RunHeightMapShape(&hm_geometry, tf1, hm1_displacement, &shape, tf2,
                             request, result);
    const auto& contacts = result.raw_contacts();
    std::set<int> contact_set;
    for (const auto& contact : contacts) {
      contact_set.insert(contact.b2);
    }

    // std::cerr << "The # of contacts in translational heightmap test "
    //           << contact_set.size() << std::endl;

    const heightmap::FlatHeightMap<S>& bottom_map =
        hm_geometry.raw_heightmap()->bottom();
    for (uint16_t y = 0; y < bottom_map.full_shape_y(); y++) {
      for (uint16_t x = 0; x < bottom_map.full_shape_x(); x++) {
        Box<S> box;
        Transform3<S> box_tf;
        const heightmap::Pixel pixel = {x, y};
        bool ok = bottom_map.pixelToBox(pixel, box, box_tf);
        if (!ok) continue;
        box.computeLocalAABB();
        using PrimitiveSolver =
            detail::ShapePairTranslationalCollisionSolver<S>;
        ContinuousCollisionResult<S> box_result;
        PrimitiveSolver::template RunShapePair<Box<S>, Shape>(
            &box, tf1 * box_tf, hm1_displacement, &shape, tf2, request,
            box_result);
        bool is_collision_with_bin = (box_result.num_contacts() > 0);
        const int code = heightmap::encodePixel(pixel);
        bool has_this_contact = (contact_set.find(code) != contact_set.end());
        EXPECT_EQ(has_this_contact, is_collision_with_bin);
      }
    }
  }
}

template <typename S>
void heightmapShapeTranslationalCollisionRandomTests(std::size_t test_n) {
  std::vector<fcl::HeightMapCollisionGeometry<S>> hm_geometries;
  {
    const auto points = generateRandomPointCloud(100000);
    auto heightMap =
        std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
    heightMap->updateHeightsByPointCloud3D(*points);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  // Test loop
  std::array<S, 6> extent{-1, -1, -1, 1, 1, 1};
  for (std::size_t test_i = 0; test_i < test_n; test_i++) {
    for (const auto& geometry : hm_geometries) {
      Eigen::Transform<S, 3, Eigen::Isometry> tf1;
      test::generateRandomTransform(extent, tf1);
      Eigen::Transform<S, 3, Eigen::Isometry> tf2;
      test::generateRandomTransform(extent, tf2);

      // Random displacement
      TranslationalDisplacement<S> hm_displacement;
      {
        Transform3<S> translation_pose;
        test::generateRandomTransform(extent, translation_pose);
        Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

        // Setup the axis
        hm_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
        hm_displacement.scalar_displacement =
            std::abs(translation_pose.translation()[0]);
      }

      fcl::Box<S> box(0.4, 0.4, 0.4);
      box.computeLocalAABB();
      heightmapShapeTranslationalCollisionRandomTest<S, Box<S>>(
          geometry, tf1, hm_displacement, box, tf2);

      fcl::Sphere<S> sphere(0.4);
      sphere.computeLocalAABB();
      heightmapShapeTranslationalCollisionRandomTest<S, Sphere<S>>(
          geometry, tf1, hm_displacement, sphere, tf2);

      fcl::Cylinder<S> cylinder(0.4, 0.5);
      cylinder.computeLocalAABB();
      heightmapShapeTranslationalCollisionRandomTest<S, Cylinder<S>>(
          geometry, tf1, hm_displacement, cylinder, tf2);
    }
  }
}

}  // namespace fcl

GTEST_TEST(HeightMapShapeTranslationalCollision, RandomTest) {
  fcl::heightmapShapeTranslationalCollisionRandomTests<float>(10);
  fcl::heightmapShapeTranslationalCollisionRandomTests<double>(10);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
