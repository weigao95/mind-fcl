//
// Created by Wei Gao on 2024/6/13.
//
#include <gtest/gtest.h>

#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/continuous_collision.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename S>
void sphereCapsuleTest(std::size_t test_n,
                       TimeOfCollisionRequestType request_type) {
  constexpr S s1_radius = 0.1;
  constexpr S s2_radius = 0.125;
  fcl::Sphere<S> s1(s1_radius);
  fcl::Sphere<S> s2(s2_radius);
  s1.computeLocalAABB();
  s2.computeLocalAABB();

  // Random pose
  const S extent_scalar = 0.5;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};
  for (std::size_t i = 0; i < test_n; i++) {
    Transform3<S> b1_pose, b2_pose;
    test::generateRandomTransform(extent, b1_pose);
    test::generateRandomTransform(extent, b2_pose);

    // Random displacement
    TranslationalDisplacement<S> s1_displacement;
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);

      // Setup the axis, only unit z
      s1_displacement.unit_axis_in_shape1 = Vector3<S>::UnitZ();
      s1_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    // Try with capsule
    Capsule<S> capsule(s1_radius, s1_displacement.scalar_displacement);
    capsule.computeLocalAABB();
    Transform3<S> capsule_in_sphere;
    capsule_in_sphere.setIdentity();
    capsule_in_sphere.translation()[2] =
        S(0.5) * s1_displacement.scalar_displacement;
    Transform3<S> capsule_pose = b1_pose * capsule_in_sphere;
    fcl::CollisionRequest<S> capsule_request;
    fcl::CollisionResult<S> capsule_result;
    fcl::collide(&capsule, capsule_pose, &s2, b2_pose, capsule_request,
                 capsule_result);
    const bool is_overlap_capsule = (capsule_result.numContacts() > 0);

    // Try with solver
    ContinuousCollisionRequest<S> request;
    request.request_type = request_type;

    ContinuousCollisionResult<S> result;
    ShapePairTranslationalCollisionSolver<S>::RunShapePair(
        &s1, b1_pose, s1_displacement, &s2, b2_pose, request, result);
    const bool is_overlap_solver = (result.num_contacts() > 0);

    // Should match
    EXPECT_EQ(is_overlap_solver, is_overlap_capsule);

    // Try with public interface
    ContinuousCollisionResult<S> result2;
    translational_ccd(&s1, b1_pose, s1_displacement, &s2, b2_pose, request,
                      result2);
    const bool is_overlap_interface = (result2.num_contacts() > 0);

    // Should match
    EXPECT_EQ(is_overlap_interface, is_overlap_capsule);

    // For one toc
    if (is_overlap_interface &&
        request_type == TimeOfCollisionRequestType::kOneTocSample) {
      EXPECT_TRUE(!result2.raw_contacts().empty());
      const auto toc_sample = result2.raw_contacts()[0].toc.lower_bound;
      EXPECT_NEAR(result2.raw_contacts()[0].toc.upper_bound, toc_sample, 1e-6);

      // Local transform
      Transform3<S> sphere_local_tf_toc;
      sphere_local_tf_toc.setIdentity();
      sphere_local_tf_toc.translation()[2] =
          toc_sample * s1_displacement.scalar_displacement;
      Transform3<S> sphere_pose_toc = b1_pose * sphere_local_tf_toc;

      // Into fcl
      const auto translation_diff =
          (b2_pose.translation() - sphere_pose_toc.translation()).norm();
      const bool is_overlap_toc = (translation_diff <= (s1.radius + s2.radius));
      if (!is_overlap_toc) {
        ContinuousCollisionResult<S> result3;
        translational_ccd(&s1, b1_pose, s1_displacement, &s2, b2_pose, request,
                          result3);
        auto new_toc = result3.raw_contacts()[0].toc.lower_bound;
        EXPECT_NEAR(new_toc, toc_sample, 1e-5);
      }
      EXPECT_TRUE(is_overlap_toc);
    }
  }
}

template <typename S>
void boxPairIntersectTest(std::size_t test_n) {
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
  box1.computeLocalAABB();
  box2.computeLocalAABB();

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

    // Compute the overlap by obb
    Interval<S> interval;
    const bool is_overlap_box = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);

    // Compute the overlap by shape interface
    ContinuousCollisionRequest<S> request;
    ContinuousCollisionResult<S> result;
    ShapePairTranslationalCollisionSolver<S>::RunShapePair(
        &box1, b1_pose, box1_displacement, &box2, b2_pose, request, result);
    const bool is_overlap_solver = (result.num_contacts() > 0);
    EXPECT_EQ(is_overlap_box, is_overlap_solver);
  }
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(ShapePairCCD, BoxTest) {
  fcl::detail::boxPairIntersectTest<float>(100000);
  fcl::detail::boxPairIntersectTest<double>(100000);
}

GTEST_TEST(ShapePairCCD, SphereCapsuleTest) {
  fcl::detail::sphereCapsuleTest<float>(
      100000, fcl::TimeOfCollisionRequestType::kNotRequested);
  fcl::detail::sphereCapsuleTest<double>(
      100000, fcl::TimeOfCollisionRequestType::kNotRequested);

  fcl::detail::sphereCapsuleTest<float>(
      100000, fcl::TimeOfCollisionRequestType::kOneTocSample);
  fcl::detail::sphereCapsuleTest<double>(
      100000, fcl::TimeOfCollisionRequestType::kOneTocSample);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}