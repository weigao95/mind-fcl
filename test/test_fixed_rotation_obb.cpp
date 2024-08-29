//
// Created by wei on 24-3-30.
//

#include <gtest/gtest.h>

#include "fcl/math/bv/OBB.h"
#include "fcl/math/fixed_rotation_obb_disjoint.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void fixedOrientationObbDisjointInstance() {
  AABB<S> aabb1;
  aabb1.min_ = Vector3<S>(0.3840000033378601074, -0.6079999804496765137, 0);
  aabb1.max_ = Vector3<S>(0.4160000085830688477, -0.5760000348091125488,
                          1.749000072479248047);

  AABB<S> aabb2;
  aabb2.min_ = Vector3<S>(-0.4480000138282775879, -0.6399999856948852539,
                          0.1439999938011169434);
  aabb2.max_ = Vector3<S>(-0.4320000112056732178, -0.6239999532699584961,
                          0.1599999964237213135);

  fcl::Transform3<S> obj1_pose, obj2_pose;
  {
    auto& shape1_pose = obj1_pose;
    shape1_pose.setIdentity();
    shape1_pose.translation().x() = -0.21502685546875;
    shape1_pose.translation().y() = 0.365966796875;
    shape1_pose.translation().z() = -0.002374269068241119385;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.4451195597648620605, 0.4847250878810882568, 0.7529343366622924805);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.8949251770973205566, 0.2114428281784057617, 0.3929387032985687256);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.03126472234725952148, -0.8487246036529541016, 0.5279099941253662109);
    shape1_pose.linear().matrix() = rotation_matrix;
  }

  // obj2 pose
  {
    auto& shape2_pose = obj2_pose;
    shape2_pose.setIdentity();
    shape2_pose.translation().x() = -0.2467041015625;
    shape2_pose.translation().y() = -0.44219970703125;
    shape2_pose.translation().z() =  0.023712158203125;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.1737782508134841919, -0.03557235002517700195, 0.9841421246528625488);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.9837378859519958496, 0.03979833424091339111, 0.1751454174518585205);
    rotation_matrix.row(2) =
        fcl::Vector3<S>(-0.04539754986763000488, -0.998574376106262207,
                        -0.02807778120040893555);
    shape2_pose.linear().matrix() = rotation_matrix;
  }

  OBB<S> obb1, obb2;
  convertBV(aabb1, obj1_pose, obb1);
  convertBV(aabb2, obj2_pose, obb2);
  bool is_overlap_1 = obb1.overlap(obb2);

  FixedRotationBoxDisjoint<S> obb_disjoint;
  obb_disjoint.initialize(obj1_pose, obj2_pose);
  bool is_disjoint_2 = obb_disjoint.isDisjoint(aabb1, aabb2, true);
  EXPECT_EQ(is_overlap_1, !is_disjoint_2);
}

template <typename S>
void fixedOrientationObbDisjointTest(std::size_t test_n) {
  fcl::AABB<S> aabb1;
  aabb1.min_.setConstant(-0.1);
  aabb1.max_.setConstant(0.2);
  fcl::AABB<S> aabb2 = aabb1;

  // Random pose
  const S extent_scalar = 0.5;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};
  std::size_t non_strict_mismatch_count = 0;
  for (std::size_t i = 0; i < test_n; i++) {
    Transform3<S> aabb1_pose, aabb2_pose;
    test::generateRandomTransform(extent, aabb1_pose);
    test::generateRandomTransform(extent, aabb2_pose);

    // Try with obb
    OBB<S> obb1, obb2;
    convertBV(aabb1, aabb1_pose, obb1);
    convertBV(aabb2, aabb2_pose, obb2);
    const bool is_overlap_1 = obb1.overlap(obb2);

    // Try with raw disjoint
    FixedRotationBoxDisjoint<S> obb_disjoint;
    obb_disjoint.initialize(aabb1_pose, aabb2_pose);
    const bool is_overlap_2 = !obb_disjoint.isDisjoint(aabb1, aabb2, true);
    EXPECT_EQ(is_overlap_1, is_overlap_2);

    const bool is_overlap_3 = !obb_disjoint.isDisjoint(aabb1, aabb2, false);
    if (is_overlap_1 != is_overlap_3) {
      non_strict_mismatch_count++;
    }
  }

  std::cout << "Mismatch count " << non_strict_mismatch_count << " among "
            << test_n << " attempts. Mismatch ratio is "
            << double(non_strict_mismatch_count) / double(test_n) << std::endl;
}

}  // namespace fcl

GTEST_TEST(FixedOrientationObbTest, CompareWithObb) {
  fcl::fixedOrientationObbDisjointInstance<float>();
  fcl::fixedOrientationObbDisjointTest<float>(100000);
  fcl::fixedOrientationObbDisjointTest<double>(100000);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
