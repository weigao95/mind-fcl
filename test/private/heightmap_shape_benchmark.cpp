#include <array>
#include <chrono>
#include <limits>

#include "fcl/geometry/heightmap/heightmap_collision_geometry.h"
#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

void buildHeightMap(
    fcl::heightmap::LayeredHeightMap<float>& heightmap_to_build) {
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return;

  // Forward to builder
  std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";
  auto point_cloud = fcl::test::readPointCloudXYZ(pcl_xyz_path);
  // Add a offset to z
  auto cloud_with_offset = std::make_shared<fcl::heightmap::PointCloud>();
  float z_offset = 0.0109;
  for(std::size_t i = 0; i < point_cloud->size(); i++) {
    auto point_i = point_cloud->getPoint(i);
    point_i.z() += z_offset;
    cloud_with_offset->push_back(point_i);
  }

  heightmap_to_build.resetHeights();
  heightmap_to_build.updateHeightsByPointCloud3D(*cloud_with_offset);
}

void heightmapSphereBenchmark() {
  // Make the cloud
  auto height_map =
      std::make_shared<fcl::heightmap::LayeredHeightMap<float>>(0.001f, 1024);
  buildHeightMap(*height_map);
  auto hm_geometry =
      std::make_shared<fcl::HeightMapCollisionGeometry<float>>(height_map);
  Eigen::Transform<float, 3, Eigen::Isometry> pose_cloud;
  // replicate the setup in octree, z offset is partially added in buildHM
  pose_cloud.setIdentity();
  // pose_cloud.translation().z() += 0.0109;
  fcl::CollisionObject<float> cloud_obj(hm_geometry, pose_cloud);

  auto sphere_geometry = std::make_shared<fcl::Sphere<float>>(0.1f);
  Eigen::Transform<float, 3, Eigen::Isometry> pose_sphere;
  pose_sphere.setIdentity();
  fcl::CollisionObject<float> sphere_obj(sphere_geometry, pose_sphere);

  // Try collide
  fcl::CollisionRequest<float> request;
  request.disablePenetration();
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
  std::cout << "Time is ms " << ms_time << " for test_n runs " << test_n << std::endl;
  std::cout << "The number of contact is " << result.numContacts() << std::endl;
}

int main() {
  heightmapSphereBenchmark();
  return 0;
}