//
// Created by Wei Gao on 2024/6/14.
//
#include <gtest/gtest.h>

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/narrowphase/detail/ccd/bvh_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename BV>
std::shared_ptr<const BVHModel<BV>> generateBoxBVHModel() {
  using S = typename BV::S;
  Box<S> box(0.4, 0.4, 0.4);
  return test::generateBoxBVHModel<BV>(box);
}

template <typename BV>
void shapeBVHSimpleBoxTest(std::size_t test_n) {
  using S = typename BV::S;
  Box<S> box(0.4, 0.4, 0.4);
  box.computeLocalAABB();
  auto box_bvh = generateBoxBVHModel<BV>();

  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (std::size_t test_i = 0; test_i < test_n; test_i++) {
    Transform3<S> tf1;
    test::generateRandomTransform(extent, tf1);
    Transform3<S> tf2;
    test::generateRandomTransform(extent, tf2);

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

    // As obb
    OBB<S> obb1, obb2;
    computeBV(box, tf1, obb1);
    computeBV(box, tf2, obb2);
    Interval<S> interval;
    const bool is_intersect_obb = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);

    ContinuousCollisionRequest<S> request;
    ContinuousCollisionResult<S> result1;
    TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunShapeMesh(
        &box, tf1, box1_displacement, box_bvh.get(), tf2, request, result1);
    bool is_intersect_solver1 = (result1.num_contacts() > 0);
    EXPECT_EQ(is_intersect_solver1, is_intersect_obb);

    // swap and try again
    TranslationalDisplacement<S> s2_displacement;
    s2_displacement.unit_axis_in_shape1 =
        tf2.linear().transpose() *
        (tf1.linear() * (-box1_displacement.unit_axis_in_shape1));
    s2_displacement.scalar_displacement = box1_displacement.scalar_displacement;

    ContinuousCollisionResult<S> result2;
    TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunMeshShape(
        box_bvh.get(), tf2, s2_displacement, &box, tf1, request, result2);
    bool is_intersect_solver2 = (result2.num_contacts() > 0);
    EXPECT_EQ(is_intersect_solver2, is_intersect_obb);
  }
}

template <typename BV>
void bvhPairSimpleBoxTest(std::size_t test_n) {
  using S = typename BV::S;
  Box<S> box(0.4, 0.4, 0.4);
  box.computeLocalAABB();
  auto box_bvh = generateBoxBVHModel<BV>();

  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (std::size_t test_i = 0; test_i < test_n; test_i++) {
    Transform3<S> tf1;
    test::generateRandomTransform(extent, tf1);
    Transform3<S> tf2;
    test::generateRandomTransform(extent, tf2);

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

    ContinuousCollisionRequest<S> request;
    ContinuousCollisionResult<S> result;
    TranslationalDisplacementOrientedNodeBVHSolver<BV>::RunMeshPair(
        box_bvh.get(), tf1, box1_displacement, box_bvh.get(), tf2, request,
        result);
    bool is_intersect_solver = (result.num_contacts() > 0);

    // As obb
    OBB<S> obb1, obb2;
    computeBV(box, tf1, obb1);
    computeBV(box, tf2, obb2);
    Interval<S> interval;
    const bool is_intersect_obb = !BoxPairTranslationalCCD<S>::IsDisjoint(
        obb1, box1_displacement, obb2, interval, 1e-4);
    EXPECT_EQ(is_intersect_solver, is_intersect_obb);
  }
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(BVH_CCD_Test, PrimitiveTest) {
  fcl::detail::shapeBVHSimpleBoxTest<fcl::AABB<float>>(100000);
  fcl::detail::shapeBVHSimpleBoxTest<fcl::AABB<double>>(100000);
  fcl::detail::shapeBVHSimpleBoxTest<fcl::OBB<float>>(100000);
  fcl::detail::shapeBVHSimpleBoxTest<fcl::OBB<double>>(100000);
}

GTEST_TEST(BVH_CCD_Test, PrimitiveBVH_PairTest) {
  fcl::detail::bvhPairSimpleBoxTest<fcl::AABB<float>>(100000);
  fcl::detail::bvhPairSimpleBoxTest<fcl::AABB<double>>(100000);
  fcl::detail::bvhPairSimpleBoxTest<fcl::OBB<float>>(100000);
  fcl::detail::bvhPairSimpleBoxTest<fcl::OBB<double>>(100000);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}