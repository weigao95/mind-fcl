//
// Created by wei on 2/2/22.
//
#include <array>
#include <chrono>
#include <limits>

#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

void octree2MeshContactBenchmark() {
  // Get the path
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return;

  // Get the model path
  std::string mesh_path = test_model_dir + "/complex_3.STL";
  std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";

  // Read the mesh
  std::vector<fcl::Vector3f> points;
  std::vector<fcl::MeshSimplex> triangles;
  fcl::test::loadSTLMesh(mesh_path, points, triangles);
  std::cout << "# points " << points.size() << std::endl;
  std::cout << "# triangles " << triangles.size() << std::endl;

  // Build bvh
  auto const bvhm = std::make_shared<fcl::BVHModel<fcl::OBBf>>();
  bvhm->beginModel();
  bvhm->addSubModel(points, triangles);
  bvhm->endModel();

  // Make the object
  Eigen::Transform<float, 3, Eigen::Isometry> pose_mesh3;
  pose_mesh3.setIdentity();
  bvhm->computeLocalAABB();
  fcl::CollisionObject<float> mesh3_obj(
      bvhm, pose_mesh3, fcl::CollisionObjectf::GeometryLocalAABBComputed());

  auto fcl_tree_0 = fcl::test::buildOctree2(pcl_xyz_path);
  Eigen::Transform<float, 3, Eigen::Isometry> pose_cloud;
  pose_cloud.setIdentity();
  pose_cloud.translation().z() += 0.0109;
  fcl::CollisionObject<float> cloud_obj(fcl_tree_0, pose_cloud);

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
  std::cout << "Time is ms " << ms_time << " for " << test_n << " runs. "
            << std::endl;
}

int main() {
  std::cout << "Octree2: " << std::endl;
  octree2MeshContactBenchmark();
}