//
// Created by wei on 2/2/22.
//
#include <array>
#include <chrono>
#include <limits>

#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

std::shared_ptr<fcl::Octree2CollisionGeometry<float>> buildOctree2(
    float voxel_resolution = 0.001, std::uint16_t bottom_half_shape = 1024) {
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return nullptr;

  // Forward to builder
  std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";
  return fcl::test::buildOctree2(pcl_xyz_path, voxel_resolution,
                                 bottom_half_shape);
}

void octree2SphereBenchmark() {
  // Make the cloud
  auto fcl_tree_0 = buildOctree2();
  Eigen::Transform<float, 3, Eigen::Isometry> pose_cloud;
  pose_cloud.setIdentity();
  pose_cloud.translation().z() += 0.0109;
  fcl::CollisionObject<float> cloud_obj(fcl_tree_0, pose_cloud);

  auto sphere_geometry = std::make_shared<fcl::Sphere<float>>(0.1f);
  Eigen::Transform<float, 3, Eigen::Isometry> pose_sphere;
  pose_sphere.setIdentity();
  fcl::CollisionObject<float> sphere_obj(sphere_geometry, pose_sphere);

  // Try collide
  fcl::CollisionRequest<float> request;
  request.setMaxContactCount(50000000);
  request.disablePenetration();
  fcl::CollisionResult<float> result;
  auto start = std::chrono::high_resolution_clock::now();
  int test_n = 400;
  {
    for (auto j = 0; j < test_n; j++) {
      result.clear();
      fcl::collide(&sphere_obj, &cloud_obj, request, result);
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "In octree2 vs shape, time is ms " << ms_time
            << " for test_n runs " << test_n << std::endl;
  std::cout << "The number of contact is " << result.numContacts() << std::endl;
}

int main() {
  octree2SphereBenchmark();
  return 0;
}