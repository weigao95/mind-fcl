//
// Created by Wei Gao on 2024/6/18.
//
#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/continuous_collision.h"
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

template <typename S>
int heightMapOctreeTranslationalCollisionCompareWithPrimitive(
    std::shared_ptr<const octree2::Octree<S>> octree,
    std::shared_ptr<const heightmap::LayeredHeightMap<S>> height_map,
    const Transform3<S>& tf_octree, const Transform3<S>& tf_hm,
    const TranslationalDisplacement<S>& hm_displacement) {
  detail::TranslationalDisplacementHeightMapSolver<S> solver;

  // Construct the geometry
  fcl::Octree2CollisionGeometry<S> octree_geo0(octree);
  octree_geo0.computeLocalAABB();
  fcl::HeightMapCollisionGeometry<S> hm_geometry(height_map);
  hm_geometry.computeLocalAABB();

  // First checking
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 1000000;
  ContinuousCollisionResult<S> result0;
  solver.RunHeightMapOctree(&hm_geometry, tf_hm, hm_displacement, &octree_geo0,
                            tf_octree, request, result0);

  ContinuousCollisionResult<S> result1;
  fcl::translational_ccd(&hm_geometry, tf_hm, hm_displacement, &octree_geo0,
                         tf_octree, request, result1);
  EXPECT_EQ(result0.num_contacts(), result1.num_contacts());

  // Check with each bv
  for (std::size_t i = 0; i < result0.num_contacts(); i++) {
    const ContinuousCollisionContact<S>& contact_i = result0.raw_contacts()[i];
    auto pixel = heightmap::decodePixel(contact_i.b1);
    AABB<S> pixel_aabb;
    auto ok = height_map->bottom().pixelToBox(pixel, pixel_aabb);
    EXPECT_TRUE(ok);

    AABB<S> octree_aabb = contact_i.o2_bv;
    OBB<S> box1, box2;
    convertBV(pixel_aabb, tf_hm, box1);
    convertBV(octree_aabb, tf_octree, box2);
    Interval<S> box_interval;
    const bool intersect_by_box =
        !detail::BoxPairTranslationalCCD<S>::IsDisjoint(
            box1, hm_displacement, box2, box_interval,
            request.zero_movement_tolerance);
    EXPECT_TRUE(intersect_by_box);
  }

  // Do naive checking
  TranslationalDisplacement<S> octree_equiv_disp;
  octree_equiv_disp.scalar_displacement = hm_displacement.scalar_displacement;
  octree_equiv_disp.unit_axis_in_shape1 =
      tf_octree.linear().transpose() * tf_hm.linear() *
      (-hm_displacement.unit_axis_in_shape1);

  std::size_t naive_count = 0;
  auto visit_octree_leaf_node = [&](const AABB<S>& aabb) -> bool {
    Box<S> box_aabb;
    Transform3<S> box_tf;
    constructBox(aabb, tf_octree, box_aabb, box_tf);
    box_aabb.computeLocalAABB();
    ContinuousCollisionResult<S> box_result;
    solver.RunShapeHeightMap(&box_aabb, box_tf, octree_equiv_disp, &hm_geometry,
                             tf_hm, request, box_result);
    naive_count += box_result.num_contacts();
    return false;
  };

  octree_geo0.visitLeafNodes(visit_octree_leaf_node);

  // Return the contact size
  EXPECT_EQ(result0.num_contacts(), naive_count);
  return result0.num_contacts();
}

template <typename S>
void heightMapOctree2PrimitiveTranslationalTest(int test_n) {
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

  // Loop over pairs
  std::vector<int> n_contact_for_each_pair;
  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (int test_idx = 0; test_idx < test_n; test_idx++) {
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

    // Record the contacts
    n_contact_for_each_pair.clear();
    for (std::size_t i = 0; i < octree_hm_pairs.size(); i++) {
      const auto& test_pair_i = octree_hm_pairs[i];
      auto octree = test_pair_i.first;
      auto height_map = test_pair_i.second;
      int n_contact =
          heightMapOctreeTranslationalCollisionCompareWithPrimitive<S>(
              octree, height_map, tf1, tf2, hm_displacement);
      n_contact_for_each_pair.push_back(n_contact);
    }
  }
}

}  // namespace fcl

GTEST_TEST(HeightMapOctreeTranslationalTest, PrimitiveTest) {
  fcl::heightMapOctree2PrimitiveTranslationalTest<float>(10);
  fcl::heightMapOctree2PrimitiveTranslationalTest<double>(10);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}