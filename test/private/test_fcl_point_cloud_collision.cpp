//
// Created by mech-mind_gw on 2021/12/28.
//

#include <gtest/gtest.h>

#include <array>

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

void octreePairTest_0() {
  using S = float;
  S resolution = 0.001;
  S delta_z = 0.0109;
  int max_n_contacts = 500000;

  {
    // Should be ok as it is const ptr inside?
    auto pcl_tree = buildOctree2(resolution);
    if (pcl_tree == nullptr) return;
    pcl_tree->computeLocalAABB();
    const auto& raw_tree = pcl_tree->raw_octree();

    // Into two trees
    std::shared_ptr<const fcl::Octree2CollisionGeometry<S>> fcl_tree_0(
        pcl_tree);
    std::shared_ptr<const fcl::Octree2CollisionGeometry<S>> fcl_tree_1(
        pcl_tree);

    // Make the object
    Eigen::Transform<S, 3, Eigen::Isometry> pose_0;
    Eigen::Transform<S, 3, Eigen::Isometry> pose_1;
    pose_0.setIdentity();
    pose_1.setIdentity();
    pose_0.translation().z() += delta_z;  // 2 [cm]
    using GeometryLocalAABBComputed =
        typename fcl::CollisionObject<S>::GeometryLocalAABBComputed;
    fcl::CollisionObject<S> tree_obj_0(fcl_tree_0, pose_0,
                                       GeometryLocalAABBComputed());
    fcl::CollisionObject<S> tree_obj_1(fcl_tree_1, pose_1,
                                       GeometryLocalAABBComputed());

    // Try collide
    fcl::CollisionRequest<S> request;
    request.setMaxContactCount(max_n_contacts);
    request.useDefaultPenetration();
    fcl::CollisionResult<S> result;
    fcl::collide(&tree_obj_0, &tree_obj_1, request, result);
    std::cerr << "The #of contacts is " << result.numContacts() << std::endl;

    // Check that they all collide
    for (std::size_t i = 0; i < result.numContacts(); i++) {
      const auto& o1_bv = result.getContact(i).o1_bv;
      const auto& o2_bv = result.getContact(i).o2_bv;
      fcl::OBB<S> obb_1, obb_2;
      fcl::convertBV(o1_bv, pose_0, obb_1);
      fcl::convertBV(o2_bv, pose_1, obb_2);
      EXPECT_TRUE(obb_1.overlap(obb_2));

      // Check we can find the in original tree
      const auto& o1_center = o1_bv.center();
      EXPECT_TRUE(raw_tree->isPointOccupied(o1_center));
      const auto& o2_center = o2_bv.center();
      EXPECT_TRUE(raw_tree->isPointOccupied(o2_center));
    }
  }
}

GTEST_TEST(OctomapPointCloudTest, BasicTest) { octreePairTest_0(); }

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}