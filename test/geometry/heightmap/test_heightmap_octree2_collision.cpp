//
// Created by Wei Gao on 2024/3/27.
//

#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/collision.h"
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

template <typename S>
int heightMapOctree2CollisionCompareWithPrimitive(
    std::shared_ptr<const octree2::Octree<S>> octree,
    std::shared_ptr<const heightmap::LayeredHeightMap<S>> height_map,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_hm) {
  detail::GJKSolver<S> narrowphase_solver;
  detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);

  // Construct the geometry
  fcl::Octree2CollisionGeometry<S> octree_geo0(octree);
  octree_geo0.computeLocalAABB();
  fcl::HeightMapCollisionGeometry<S> hm_geometry(height_map);
  hm_geometry.computeLocalAABB();

  // First checking
  CollisionRequest<S> request{UINT_MAX};
  CollisionResult<S> result0;
  solver.HeightMapOctreeIntersect(&hm_geometry, &octree_geo0, tf_hm, tf_octree,
                                  request, result0);

  CollisionResult<S> result1;
  fcl::collide(&octree_geo0, tf_octree, &hm_geometry, tf_hm, request, result1);
  EXPECT_EQ(result0.numContacts(), result1.numContacts());

  // Check with each bv
  for (std::size_t i = 0; i < result0.numContacts(); i++) {
    const Contact<S>& contact_i = result0.getContact(i);
    auto pixel = heightmap::decodePixel(contact_i.b1);
    AABB<S> pixel_aabb;
    auto ok = height_map->bottom().pixelToBox(pixel, pixel_aabb);
    EXPECT_TRUE(ok);

    AABB<S> octree_aabb = contact_i.o2_bv;
    Box<S> box1, box2;
    Transform3<S> box1_tf, box2_tf;
    constructBox(pixel_aabb, tf_hm, box1, box1_tf);
    constructBox(octree_aabb, tf_octree, box2, box2_tf);
    bool intersect_by_primitive = narrowphase_solver.shapeIntersect(
        box1, box1_tf, box2, box2_tf, nullptr);
    EXPECT_TRUE(intersect_by_primitive);
  }

  // Do naive checking
  std::size_t naive_count = 0;
  auto visit_octree_leaf_node = [&](const AABB<S>& aabb) -> bool {
    Box<S> box_aabb;
    Transform3<S> box_tf;
    constructBox(aabb, tf_octree, box_aabb, box_tf);
    result1.clear();
    solver.ShapeHeightMapIntersect(box_aabb, &hm_geometry, box_tf, tf_hm,
                                   request, result1);
    naive_count += result1.numContacts();
    return false;
  };

  octree_geo0.visitLeafNodes(visit_octree_leaf_node);

  // Return the contact size
  EXPECT_EQ(result0.numContacts(), naive_count);
  return result0.numContacts();
}

template <typename S>
void heightMapOctree2PrimitiveTest() {
  using OctreeHeightMapPair =
      std::pair<std::shared_ptr<const octree2::Octree<S>>,
                std::shared_ptr<const heightmap::LayeredHeightMap<S>>>;
  std::vector<OctreeHeightMapPair> octree_hm_pairs;

  auto point_cloud_0 = generateRandomPointCloud(10000);
  auto point_fn = [point_cloud_0](int index, S& x, S& y, S& z) -> void {
    const auto point = point_cloud_0->getPoint(index);
    x = point.x();
    y = point.y();
    z = point.z();
  };

  {
    S voxel_resolution = 0.002;
    Vector3<S> resolution(voxel_resolution, voxel_resolution, voxel_resolution);
    auto octree0 = std::make_shared<octree2::Octree<S>>(resolution, 1024);
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
    octree0->rebuildTree(point_fn, point_cloud_0->size());
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
  }

  {
    S voxel_resolution = 0.002;
    Vector3<S> resolution(voxel_resolution, voxel_resolution, voxel_resolution);
    auto octree0 = std::make_shared<octree2::Octree<S>>(resolution, 1024);
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.004, 256);
    octree0->rebuildTree(point_fn, point_cloud_0->size());
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
  }

  {
    S voxel_resolution = 0.002;
    Vector3<S> resolution(voxel_resolution, voxel_resolution, voxel_resolution);
    auto octree0 = std::make_shared<octree2::Octree<S>>(resolution, 1024);
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 0.004,
                                                               512, 256);
    octree0->rebuildTree(point_fn, point_cloud_0->size());
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
    octree_hm_pairs.emplace_back(std::make_pair(octree0, hm));
  }

  // Loop over pairs
  int n_loop_per_pair = 10;
  std::vector<int> n_contact_for_each_pair;
  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (int test_idx = 0; test_idx < n_loop_per_pair; test_idx++) {
    Eigen::Transform<S, 3, Eigen::Isometry> tf1;
    test::generateRandomTransform(extent, tf1);
    Eigen::Transform<S, 3, Eigen::Isometry> tf2;
    test::generateRandomTransform(extent, tf2);

    // Record the contacts
    n_contact_for_each_pair.clear();
    for (std::size_t i = 0; i < octree_hm_pairs.size(); i++) {
      const auto& test_pair_i = octree_hm_pairs[i];
      auto octree = test_pair_i.first;
      auto height_map = test_pair_i.second;
      int n_contact = heightMapOctree2CollisionCompareWithPrimitive<S>(
          octree, height_map, tf1, tf2);
      n_contact_for_each_pair.push_back(n_contact);
    }
  }
}

}  // namespace fcl

GTEST_TEST(HeightMapOctree2Test, PrimitiveTest) {
  // std::cout << "Test with float" << std::endl;
  fcl::heightMapOctree2PrimitiveTest<float>();
  // std::cout << "Test with double" << std::endl;
  fcl::heightMapOctree2PrimitiveTest<double>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}