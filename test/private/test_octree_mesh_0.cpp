//
// Created by wei on 2/2/22.
//
#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <limits>

#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

int octreeMeshContactCount() {
  // return 3093;
  // Get the path
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model)
    return -1;

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

  // Make the cloud
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
  fcl::CollisionResult<float> local_result;
  fcl::collide(&mesh3_obj, &cloud_obj, request, local_result);
  return local_result.numContacts();
}

GTEST_TEST(OctreeMeshCollision, Test_0) {
  auto contact_0 = octreeMeshContactCount();

  // The # of contacts change as we use new GJK/EPA/MPR
  // constexpr int reference_contact_number = 2655;
  // for mac which use different Eigen
  constexpr int reference_contact_number = 3093;
  std::cerr << "Contact count " << contact_0 << std::endl;
  EXPECT_TRUE(std::abs(contact_0 - reference_contact_number) < 50);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}