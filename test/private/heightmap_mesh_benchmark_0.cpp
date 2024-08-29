//
// Created by mech-mind_gw on 2/23/2023.
//
#include <array>
#include <chrono>
#include <limits>

#include "fcl/geometry/heightmap/heightmap_collision_geometry.h"
#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

void buildHeightMapWithOffset(
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

void heightMapMeshContactBenchmark() {
  // Get the path
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return;

  // Get the model path
  std::string mesh_path = test_model_dir + "/complex_3.STL";

  auto const bvhm = std::make_shared<fcl::BVHModel<fcl::OBBf>>();
  {
    // Read the mesh
    std::vector<fcl::Vector3f> points;
    std::vector<fcl::MeshSimplex> triangles;
    fcl::test::loadSTLMesh(mesh_path, points, triangles);
    std::cout << "# points " << points.size() << std::endl;
    std::cout << "# triangles " << triangles.size() << std::endl;

    bvhm->beginModel();
    bvhm->addSubModel(points, triangles);
    bvhm->endModel();
  }

  // Make the object
  Eigen::Transform<float, 3, Eigen::Isometry> pose_mesh3;
  pose_mesh3.setIdentity();
  bvhm->computeLocalAABB();
  fcl::CollisionObject<float> mesh3_obj(
      bvhm, pose_mesh3, fcl::CollisionObjectf::GeometryLocalAABBComputed());

  // Make the cloud
  auto height_map =
      std::make_shared<fcl::heightmap::LayeredHeightMap<float>>(0.002, 512);
  buildHeightMapWithOffset(*height_map);
  auto hm_geometry =
      std::make_shared<fcl::HeightMapCollisionGeometry<float>>(height_map);
  Eigen::Transform<float, 3, Eigen::Isometry> pose_cloud;
  pose_cloud.setIdentity();
  fcl::CollisionObject<float> cloud_obj(hm_geometry, pose_cloud);

  // Try collide
  fcl::CollisionRequest<float> request;
  request.setMaxContactCount(50000000);
  request.disablePenetration();

  // Run detection
  int test_n = 100;
  auto start = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < test_n; i++) {
    fcl::CollisionResult<float> local_result;
    fcl::collide(&mesh3_obj, &cloud_obj, request, local_result);
    if (i == 0) {
      std::cout << "# of contacts " << local_result.numContacts() << std::endl;
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Time is ms " << ms_time / test_n << std::endl;
}

int main() {
  heightMapMeshContactBenchmark();
}