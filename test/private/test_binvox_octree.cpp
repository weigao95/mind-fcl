//
// Created by wei on 2/2/22.
//
#include <array>
#include <chrono>
#include <limits>
#include <gtest/gtest.h>

#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

void test_double_tree() {
  //Get the path
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model)
    return;

  // Get the model path
  std::string binvox_path = test_model_dir + "/81D.binvox";
  std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";

  // Build the cloud/object octree
  float resolution = 0.001;
  float delta_z = 0.0;
  std::uint16_t bottom_half_shape = 1024;
  auto tree_cloud = fcl::test::buildOctree2(pcl_xyz_path, resolution,
                                            bottom_half_shape);
  if (tree_cloud == nullptr)
    return;

  // Load tree object
  auto raw_tree_object =
      fcl::test::loadBinvoxAsOctree2(binvox_path, 1.0, bottom_half_shape);
  auto tree_object = std::make_shared<fcl::Octree2CollisionGeometry<float>>(raw_tree_object);
  if (tree_object == nullptr)
    return;

  {
    // Make fcl
    auto fcl_tree_0 = tree_cloud;
    auto fcl_tree_1 = tree_object;
    fcl_tree_0->computeLocalAABB();
    fcl_tree_1->computeLocalAABB();

    // Make the object
    Eigen::Isometry3f pose_0;
    Eigen::Isometry3f pose_1;
    pose_0.setIdentity();
    pose_1.setIdentity();
    pose_0.translation().z() += delta_z;  // 2 [cm]
    fcl::CollisionObjectf tree_obj_0(fcl_tree_0, pose_0);
    fcl::CollisionObjectf tree_obj_1(fcl_tree_1, pose_1);

    // Try collide
    fcl::CollisionRequestf request;
    request.setMaxContactCount(500000);
    request.useDefaultPenetration();
    fcl::CollisionResultf local_result;
    fcl::collide(&tree_obj_0, &tree_obj_1, request, local_result);
    const std::size_t reference_contact_count = 86596;
    EXPECT_EQ(local_result.numContacts(), reference_contact_count);
  }
}

GTEST_TEST(BinvoxTest, DoubleTree) {
  test_double_tree();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}