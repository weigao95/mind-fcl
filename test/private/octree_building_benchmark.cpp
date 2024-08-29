//
// Created by mech-mind_gw on 2/10/2022.
//

#include <array>
#include <chrono>
#include <limits>

#include "fcl/geometry/octree2/octree.h"
#include "test_fcl_utility.h"

void buildOctreeBenchmark(
    const std::shared_ptr<const fcl::octomap::Pointcloud>& point_cloud,
    double voxel_resolution = 0.001) {
  // Point fn
  auto point_fn = [&point_cloud](int index, float& x, float& y,
                                 float& z) -> void {
    const auto point = point_cloud->getPoint(index);
    x = point.x();
    y = point.y();
    z = point.z();
  };

  {
    auto start = std::chrono::high_resolution_clock::now();
    auto tree2 =
        std::make_shared<fcl::octree2::Octree<float>>(voxel_resolution, 2048);
    tree2->rebuildTree(point_fn, point_cloud->size());
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "Build octree2 time is ms " << ms_time << std::endl;

    auto tree_geom =
        std::make_shared<fcl::Octree2CollisionGeometry<float>>(tree2);
    fcl::AABB<float> prune_aabb;
    prune_aabb.max_ = fcl::Vector3<float>(0.2, 0.2, 0.2);
    prune_aabb.min_ = -prune_aabb.max_;

    // Random pose
    const float extent_scalar = 0.1;
    std::array<float, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                                extent_scalar,  extent_scalar,  extent_scalar};
    fcl::Transform3f obb_pose;
    fcl::test::generateRandomTransform(extent, obb_pose);
    fcl::OBB<float> prune_obb;
    fcl::convertBV(prune_aabb, obb_pose, prune_obb);

    // Test prune
    start = std::chrono::high_resolution_clock::now();
    auto new_tree_geom = tree_geom->pruneBy(prune_obb, true);
    end = std::chrono::high_resolution_clock::now();
    ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "Prune octree2 and rebuild time is ms " << ms_time
              << std::endl;

    // Test prune
    start = std::chrono::high_resolution_clock::now();
    new_tree_geom = tree_geom->pruneBy(prune_obb, false);
    end = std::chrono::high_resolution_clock::now();
    ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "Prune octree2 time is ms " << ms_time << std::endl;
  }
}

void octreeBuildingBenchmark() {
  // Test large cloud 0
  const double voxel_resolution = 0.001;
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) throw std::runtime_error("Cannot find the model path");

  // Default cloud
  {
    std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";
    auto large_cloud_0 = fcl::test::readPointCloudXYZ(pcl_xyz_path);
    std::cout << "A medium (in dimension) point cloud with # of points "
              << large_cloud_0->size() << std::endl;
    buildOctreeBenchmark(large_cloud_0, voxel_resolution);
  }

  // Large test 0
  {
    std::string pcl_xyz_path = test_model_dir + "/pcl_large_test_0.xyz";
    auto large_cloud_0 = fcl::test::readPointCloudXYZ(pcl_xyz_path);
    std::cout << "A small (in dimension) point cloud with # of points "
              << large_cloud_0->size() << std::endl;
    buildOctreeBenchmark(large_cloud_0, voxel_resolution);
  }

  // Large test 1
  {
    std::string pcl_xyz_path = test_model_dir + "/pcl_large_test_1.xyz";
    auto large_cloud_0 = fcl::test::readPointCloudXYZ(pcl_xyz_path);
    std::cout << "A large (in dimension) point cloud with # of points "
              << large_cloud_0->size() << std::endl;
    buildOctreeBenchmark(large_cloud_0, 0.0015);
  }
}

int main() {
  octreeBuildingBenchmark();
  return 0;
}