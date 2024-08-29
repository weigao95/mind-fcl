//
// Created by mech-mind_gw on 2/22/2023.
//

#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace heightmap {

template <typename S>
void updateByRandomPointCloud(LayeredHeightMap<S>& heightmap) {
  PointCloud point_cloud;
  std::size_t n_points = 10000;
  point_cloud.clear();
  point_cloud.reserve(n_points);
  bool exclude_out_range_points = true;
  while (point_cloud.size() < n_points) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    // Remove the point if it is too far
    if (exclude_out_range_points &&
        (std::abs(point.x()) >= heightmap.bottom().half_range_x() ||
         std::abs(point.y()) >= heightmap.bottom().half_range_y())) {
      continue;
    }

    // Some of the points are below check xOy plane
    const S z_height = point[2] + 0.9;
    point_cloud.push_back(point[0], point[1], z_height);
  }

  // Perform update
  heightmap.updateHeightsByPointCloud3D(point_cloud);
}

template <typename S>
void heightmapPairCollisionCompareWithNaive(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  // Construct the solver
  const fcl::detail::GJKSolver<S> narrowphase_solver;
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);

  const fcl::CollisionRequest<S> request{UINT_MAX};
  fcl::CollisionResult<S> result;
  solver.HeightMapIntersect(&hm_geometry, &hm_geometry, tf1, tf2, request,
                            result);
  const auto& contacts = result.getContacts();
  // std::cerr << "# contacts in heightmap pair " << contacts.size() <<
  // std::endl;
  std::set<std::pair<int, int>> contact_set;
  for (const auto& contact : contacts) {
    contact_set.insert(std::make_pair(contact.b1, contact.b2));
  }

  // Compare with fcl
  {
    CollisionResult<S> fcl_collide_result;
    fcl::collide(&hm_geometry, tf1, &hm_geometry, tf2, request,
                 fcl_collide_result);
    EXPECT_EQ(fcl_collide_result.numContacts(), result.numContacts());
  }

  // The same map
  const FlatHeightMap<S>& hm1 = hm_geometry.raw_heightmap()->bottom();
  const FlatHeightMap<S>& hm2 = hm_geometry.raw_heightmap()->bottom();
  for (uint16_t hm1_y = 0; hm1_y < hm1.full_shape_y(); hm1_y++) {
    for (uint16_t hm1_x = 0; hm1_x < hm1.full_shape_x(); hm1_x++) {
      for (uint16_t hm2_y = 0; hm2_y < hm2.full_shape_y(); hm2_y++) {
        for (uint16_t hm2_x = 0; hm2_x < hm2.full_shape_x(); hm2_x++) {
          Pixel hm1_pixel(hm1_x, hm1_y);
          Pixel hm2_pixel(hm2_x, hm2_y);

          Box<S> box1, box2;
          Transform3<S> box_tf_1, box_tf_2;
          bool has_box1 = hm1.pixelToBox(hm1_pixel, box1, box_tf_1);
          bool has_box2 = hm2.pixelToBox(hm2_pixel, box2, box_tf_2);
          box_tf_1 = tf1 * box_tf_1;
          box_tf_2 = tf2 * box_tf_2;
          bool is_collision_with_bin =
              has_box1 && has_box2 &&
              narrowphase_solver.shapeIntersect(box1, box_tf_1, box2, box_tf_2,
                                                nullptr);
          int code1 = heightmap::encodePixel(hm1_pixel);
          int code2 = heightmap::encodePixel(hm2_pixel);
          bool has_this_contact = (contact_set.find(std::make_pair(
                                       code1, code2)) != contact_set.end());
          EXPECT_EQ(has_this_contact, is_collision_with_bin);
        }
      }
    }
  }
}

template <typename S>
void heightmapPairCollisionCompareWithNaiveTest() {
  std::vector<fcl::HeightMapCollisionGeometry<S>> hm_geometries;

  // Roughly 1m scale
  {
    // No geometry
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.12, 8);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  auto update_heightmap_visitor =
      [](const Pixel& pixel, const Point2D<S>& box_bottom_center,
         uint16_t old_height_in_mm, uint16_t& new_height_in_mm) -> bool {
    (void)(box_bottom_center);
    (void)(old_height_in_mm);
    new_height_in_mm = 10 * (pixel.x + pixel.y);
    return false;
  };

  {
    // x + y
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.12, 8);
    height_map->updateHeightsByBottomLayerUpdateFunctor(
        update_heightmap_visitor);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // x + y
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.6, 16);
    height_map->updateHeightsByBottomLayerUpdateFunctor(
        update_heightmap_visitor);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // x + y
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.6, 0.12, 16, 8);
    height_map->updateHeightsByBottomLayerUpdateFunctor(
        update_heightmap_visitor);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // Random point cloud
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.12, 8);
    updateByRandomPointCloud(*height_map);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // Random point cloud
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.6, 16);
    updateByRandomPointCloud(*height_map);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // Random point cloud
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.3, 32);
    updateByRandomPointCloud(*height_map);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    // Random point cloud
    auto height_map = std::make_shared<LayeredHeightMap<S>>(0.3, 0.4, 32, 16);
    updateByRandomPointCloud(*height_map);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(height_map);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  // Test loop
  std::array<S, 6> extent{-1, -1, -0.2, 1, 1, 0.2};
  for (const auto& geometry : hm_geometries) {
    Eigen::Transform<S, 3, Eigen::Isometry> tf1;
    test::generateRandomTransform(extent, tf1);
    Eigen::Transform<S, 3, Eigen::Isometry> tf2;
    test::generateRandomTransform(extent, tf2);
    heightmapPairCollisionCompareWithNaive<S>(geometry, tf1, tf2);
  }
}

}  // namespace heightmap
}  // namespace fcl

GTEST_TEST(HeightMapPairCollisionTest, CompareWithNaive) {
  fcl::heightmap::heightmapPairCollisionCompareWithNaiveTest<float>();
  fcl::heightmap::heightmapPairCollisionCompareWithNaiveTest<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
