//
// Created by mech-mind_gw on 3/22/2024.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_util.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace octree2 {

void octreeNodeBitsetTest() {
  Bitset8 bitset;  // Default construct as empty
  EXPECT_TRUE(bitset.is_all_cleared());

  for (std::uint8_t i = 0; i < 8; i++) {
    EXPECT_FALSE(bitset.test_i(i));
  }

  for (std::uint8_t i = 0; i < 8; i++) {
    EXPECT_TRUE(bitset.count_of_bits() == i);
    bitset.set_i(i);
    EXPECT_TRUE(bitset.test_i(i));
    EXPECT_TRUE(bitset.count_of_bits() == i + 1);
  }

  EXPECT_TRUE(bitset.is_all_set());

  for (std::uint8_t i = 0; i < 8; i++) {
    EXPECT_TRUE(bitset.count_of_bits() == 8 - i);
    bitset.clear_i(i);
    EXPECT_FALSE(bitset.test_i(i));
    EXPECT_TRUE(bitset.count_of_bits() == 7 - i);
  }
}

template <typename S>
void octreeNodeContainmentTest(std::size_t test_n) {
  const S aabb_box_size = 0.1;
  const S obb_box_size = 0.3;
  AABB<S> obb_from_aabb;
  obb_from_aabb.min_.setConstant(-obb_box_size);
  obb_from_aabb.max_.setConstant(obb_box_size);

  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  Transform3<S> obb_pose;
  for (std::size_t i = 0; i < test_n; i++) {
    Vector3<S> aabb_center_i;
    aabb_center_i.setRandom();
    aabb_center_i *= 0.299;
    /*aabb_center_i.x() = -0.2653833627700805664;
    aabb_center_i.y() = -0.2937439680099487305;
    aabb_center_i.z() = 0.2504365444183349609;*/

    AABB<S> aabb_i;
    aabb_i.min_ = aabb_center_i;
    aabb_i.max_ = aabb_center_i;
    for (auto k = 0; k < 3; k++) {
      aabb_i.min_[k] -= S(0.5) * aabb_box_size;
      aabb_i.max_[k] += S(0.5) * aabb_box_size;
    }

    OBB<S> obb_i;
    test::generateRandomTransform(extent, obb_pose);

    // obj2 pose
    /*{
      auto& shape2_pose = obb_pose;
      shape2_pose.setIdentity();
      shape2_pose.translation().x() = -0.2946533262729644775;
      shape2_pose.translation().y() = -0.07327881455421447754;
      shape2_pose.translation().z() = 0.01898804306983947754;

      // The rotation matrix
      fcl::Matrix3<S> rotation_matrix;
      rotation_matrix.row(0) =
          fcl::Vector3<S>(0.7236182093620300293, -0.3470123708248138428,
                          -0.5966230630874633789);
      rotation_matrix.row(1) =
          fcl::Vector3<S>(0.002632170915603637695, 0.8658011555671691895,
                          -0.5003812909126281738);
      rotation_matrix.row(2) = fcl::Vector3<S>(
          0.6901954412460327148, 0.3605146110057830811, 0.6274228692054748535);
      shape2_pose.linear().matrix() = rotation_matrix;
    }*/

    // Test containment
    fcl::convertBV(obb_from_aabb, obb_pose, obb_i);
    bool test_contained_naive = is_contained_naive(obb_i, aabb_i);
    bool test_contained = is_contained(obb_i, aabb_i);
    EXPECT_EQ(test_contained_naive, test_contained);
  }
}

}  // namespace octree2
}  // namespace fcl

GTEST_TEST(Octree2_NodeTest, BitsetTest) {
  fcl::octree2::octreeNodeBitsetTest();
}

GTEST_TEST(Octree_NodeTest, ContainmentTest) {
  fcl::octree2::octreeNodeContainmentTest<float>(10000);
  fcl::octree2::octreeNodeContainmentTest<double>(10000);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
