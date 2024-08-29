//
// Created by Wei Gao on 2024/6/11.
//
#include <gtest/gtest.h>

#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/detail/ccd/box_pair_ccd.h"
#include "fcl/narrowphase/detail/ccd/gjk_ccd.h"
#include "fcl/narrowphase/detail/gjk_solver_cvx.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename S>
void testTocIntervalOnAxisRandom() {
  // Static case
  const S box1_half_size{0.5};
  const S box2_half_size{0.8};

  const int test_n = 100000;
  for (auto i = 0; i < test_n; i++) {
    Eigen::Vector2<S> random_vec2;
    random_vec2.setRandom();
    random_vec2 *= 10.0;

    // Compute explicit case
    const S box1_displacement = random_vec2.x();
    const S box2_offset = random_vec2.y();

    // Case 1
    Interval<S> interval_1;
    bool is_empty_1 = BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
        box1_half_size, box1_displacement, box2_half_size, box2_offset,
        interval_1);

    // Should match
    if (is_empty_1) {
      // Emtpy, check lower bound and upper bound
      const auto box2_lb = box2_offset - box2_half_size;
      const auto box2_ub = box2_offset + box2_half_size;

      // Check interval lower bound
      auto box1_lb = -box1_half_size;
      auto box1_ub = box1_half_size;
      EXPECT_TRUE(box1_lb > box2_ub || box2_lb > box1_ub);

      // Check interval upper bound
      box1_lb = box1_displacement - box1_half_size;
      box1_ub = box1_displacement + box1_half_size;
      EXPECT_TRUE(box1_lb > box2_ub || box2_lb > box1_ub);
    } else {
      // Not emtpy, check lower bound and upper bound
      const auto box2_lb = box2_offset - box2_half_size;
      const auto box2_ub = box2_offset + box2_half_size;
      const auto sign_disp = box1_displacement > 0 ? 1.0 : -1.0;

      // Check interval lower bound
      auto box1_lb = sign_disp * interval_1.lower_bound - box1_half_size;
      auto box1_ub = sign_disp * interval_1.lower_bound + box1_half_size;
      EXPECT_FALSE(box1_lb > (box2_ub + 1e-4) || box2_lb > (box1_ub + 1e-4));

      // Check interval upper bound
      box1_lb = sign_disp * interval_1.upper_bound - box1_half_size;
      box1_ub = sign_disp * interval_1.upper_bound + box1_half_size;
      EXPECT_FALSE(box1_lb > (box2_ub + 1e-4) || box2_lb > (box1_ub + 1e-4));
    }
  }
}

template <typename S>
void testTocIntervalOnAxisSymmetry() {
  // Static case
  const S box1_half_size{0.5};
  const S box2_half_size{0.8};

  // Center case with box2 offset == 0
  const int test_n = 100000;
  for (auto i = 0; i < test_n; i++) {
    Eigen::Vector2<S> random_vec2;
    random_vec2.setRandom();
    random_vec2 *= 10.0;

    // Compute explicit case
    const S box1_displacement = random_vec2.x();
    Interval<S> interval;
    bool is_empty = BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
        box1_half_size, box1_displacement, box2_half_size, 0.0, interval);
    EXPECT_FALSE(is_empty);
    EXPECT_NEAR(interval.lower_bound, 0, 1e-6);
    const auto explicit_interval_ub =
        std::min(std::abs(box1_displacement), box1_half_size + box2_half_size);
    EXPECT_NEAR(interval.upper_bound, explicit_interval_ub, 1e-6);
  }

  // Symmetry
  for (auto i = 0; i < test_n; i++) {
    Eigen::Vector2<S> random_vec2;
    random_vec2.setRandom();
    random_vec2 *= 10.0;

    // Compute explicit case
    const S box1_displacement = random_vec2.x();
    const S box2_offset = random_vec2.y();

    // Case 1
    Interval<S> interval_1;
    bool is_empty_1 = BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
        box1_half_size, box1_displacement, box2_half_size, box2_offset,
        interval_1);

    // Case 2
    Interval<S> interval_2;
    bool is_empty_2 = BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
        box1_half_size, -box1_displacement, box2_half_size, -box2_offset,
        interval_2);

    // Should match
    EXPECT_EQ(is_empty_1, is_empty_2);
    EXPECT_NEAR(interval_1.lower_bound, interval_2.lower_bound, 1e-6);
    EXPECT_NEAR(interval_1.upper_bound, interval_2.upper_bound, 1e-6);

    // Array case
    {
      Eigen::Array3<S> box1_half_size_array{box1_half_size};
      Eigen::Array3<S> box1_displacement_array{box1_displacement};
      Eigen::Array3<S> box2_half_size_array{box2_half_size};
      Eigen::Array3<S> box2_offset_array{box2_offset};
      Interval<S> interval_array;
      bool is_empty_interval = BoxPairTranslationalCCD<S>::CheckEmptyTocIntervalOnAxis(
          box1_half_size_array, box1_displacement_array, box2_half_size_array,
          box2_offset_array, interval_array, 1e-4);
      EXPECT_EQ(is_empty_1, is_empty_interval);
    }
  }
}

template <typename S>
void boxPairIntersectTest(std::size_t test_n) {
  fcl::AABB<S> aabb1;
  aabb1.min_.setConstant(-0.1);
  aabb1.max_.setConstant(0.2);
  fcl::AABB<S> aabb2 = aabb1;

  // Random pose
  const S extent_scalar = 0.5;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};
  for (std::size_t i = 0; i < test_n; i++) {
    Transform3<S> aabb1_pose, aabb2_pose;
    test::generateRandomTransform(extent, aabb1_pose);
    test::generateRandomTransform(extent, aabb2_pose);

    // Try with obb static test
    OBB<S> obb1, obb2;
    convertBV(aabb1, aabb1_pose, obb1);
    convertBV(aabb2, aabb2_pose, obb2);
    const bool is_overlap_1 = obb1.overlap(obb2);

    TranslationalDisplacement<S> box1_displacement;
    box1_displacement.unit_axis_in_shape1 = Vector3<S>::UnitX();
    box1_displacement.scalar_displacement = 0;
    Interval<S> interval;
    const bool is_overlap_2 = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);

    // Match when no displacement
    EXPECT_EQ(is_overlap_1, is_overlap_2);

    // Try a random axis
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);
      Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

      // Setup the axis
      box1_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
      box1_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    Vector3<S> unit_axis_in_world =
        obb1.axis * box1_displacement.unit_axis_in_shape1;
    EXPECT_NEAR(box1_displacement.unit_axis_in_shape1.norm(), 1.0, 1e-6);
    EXPECT_NEAR(unit_axis_in_world.norm(), 1.0, 1e-6);

    // Compute the overlap again
    const bool is_overlap_3 = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);
    if (is_overlap_3) {
      EXPECT_TRUE(interval.upper_bound >= interval.lower_bound);
      auto obb1_displaced = obb1;

      // Should intersect at lower
      const auto scalar_displacement_lower =
          box1_displacement.scalar_displacement * interval.lower_bound;
      obb1_displaced.To =
          obb1.To + (unit_axis_in_world * scalar_displacement_lower);
      EXPECT_TRUE(obb1_displaced.overlap(obb2));

      // Should intersect at upper
      const auto scalar_displacement_upper =
          box1_displacement.scalar_displacement * interval.upper_bound;
      obb1_displaced.To =
          obb1.To + (unit_axis_in_world * scalar_displacement_upper);
      EXPECT_TRUE(obb1_displaced.overlap(obb2));
    } else {
      // Test with discrete points
      const int n_discrete_points = 10;
      const S ratio = 1.0 / S(n_discrete_points);
      auto obb1_displaced = obb1;
      for (auto i = 0; i < n_discrete_points; i++) {
        const auto scale_i = i * ratio;
        const auto scalar_displacement_i =
            box1_displacement.scalar_displacement * scale_i;
        obb1_displaced.To =
            obb1.To + (unit_axis_in_world * scalar_displacement_i);
        EXPECT_FALSE(obb1_displaced.overlap(obb2));
      }
    }
  }
}

template <typename S>
void boxPairIntersectTestCompareMPR(std::size_t test_n) {
  constexpr S b1_size = 0.2;
  constexpr S b2_size = 0.25;
  constexpr S b1_half_size = b1_size * 0.5;
  constexpr S b2_half_size = b2_size * 0.5;
  fcl::AABB<S> aabb1;
  aabb1.min_.setConstant(-b1_half_size);
  aabb1.max_.setConstant(b1_half_size);
  fcl::AABB<S> aabb2;
  aabb2.min_.setConstant(-b2_half_size);
  aabb2.max_.setConstant(b2_half_size);

  fcl::Box<S> box1(b1_size, b1_size, b1_size);
  fcl::Box<S> box2(b2_size, b2_size, b2_size);
  auto gjk_box1 = detail::constructGJKGeometry(&box1);
  auto gjk_box2 = detail::constructGJKGeometry(&box2);

  // Random pose
  const S extent_scalar = 0.5;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};
  for (std::size_t i = 0; i < test_n; i++) {
    Transform3<S> b1_pose, b2_pose;
    test::generateRandomTransform(extent, b1_pose);
    test::generateRandomTransform(extent, b2_pose);

    // As obb
    OBB<S> obb1, obb2;
    convertBV(aabb1, b1_pose, obb1);
    convertBV(aabb2, b2_pose, obb2);

    // Random displacement
    TranslationalDisplacement<S> box1_displacement;
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);
      Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

      // Setup the axis
      box1_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
      box1_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    // Compute axis in world
    Vector3<S> unit_axis_in_world =
        b1_pose.linear() * box1_displacement.unit_axis_in_shape1;
    EXPECT_NEAR(box1_displacement.unit_axis_in_shape1.norm(), 1.0, 1e-6);
    EXPECT_NEAR(unit_axis_in_world.norm(), 1.0, 1e-6);

    // Compute the overlap again
    Interval<S> interval;
    const bool is_overlap_ccd = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);

    // Now try with mk diff
    const bool is_overlap_mpr =
        TranslationalCollisionGJK<S>::CheckSweptVolumeCollisionBinary(
            gjk_box1, b1_pose, box1_displacement, gjk_box2, b2_pose);

    // Should be the same
    EXPECT_EQ(is_overlap_ccd, is_overlap_mpr);

    // Add fixed orientation test
    FixedOrientationBoxPairTranslationalCCD<S> fixed_orientation_test;
    fixed_orientation_test.Initialize(b1_pose, b2_pose, box1_displacement);
    Interval<S> fixed_orientation_interval;
    bool is_overlap_fixed_orientation = !fixed_orientation_test.IsDisjoint(
        aabb1, aabb2, fixed_orientation_interval, 1e-4);
    EXPECT_EQ(is_overlap_fixed_orientation, is_overlap_mpr);
  }
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(BoxPairCCD, TocIntervalTest) {
  fcl::detail::testTocIntervalOnAxisSymmetry<float>();
  fcl::detail::testTocIntervalOnAxisSymmetry<double>();
  fcl::detail::testTocIntervalOnAxisRandom<float>();
  fcl::detail::testTocIntervalOnAxisRandom<double>();
}

GTEST_TEST(BoxPairCCD, BoxRandomTest) {
  fcl::detail::boxPairIntersectTest<float>(100000);
  fcl::detail::boxPairIntersectTest<double>(100000);
}

GTEST_TEST(BoxPairCCD, BoxGJKTest) {
  fcl::detail::boxPairIntersectTestCompareMPR<float>(100000);
  fcl::detail::boxPairIntersectTestCompareMPR<double>(100000);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}