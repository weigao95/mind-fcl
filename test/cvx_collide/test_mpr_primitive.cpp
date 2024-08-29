//
// Created by wei on 2/8/22.
//

#include <gtest/gtest.h>

#include <memory>

#include "fcl/cvx_collide/mpr.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

using namespace fcl;

// Test convert2SimplexToTetrahedron function.
// We construct a test scenario that two boxes are on the xy plane of frame F.
// The two boxes penetrate to each other, as shown in the bottom plot.
//              y
//          ┏━━━│━━━┓Box1
//         ┏┃━━┓    ┃
//      ───┃┃──┃O───┃─────x
//     box2┗┃━━┛│   ┃
//          ┗━━━│━━━┛
//
// @param[in] X_WF The pose of the frame F measured and expressed in the world
// frame W.
// @param[out] o1 Box1
// @param[out] o2 Box2
// @param[out] ccd The ccd solver info.
// @param[out] X_FB1 The pose of box 1 frame measured and expressed in frame F.
// @param[out] X_FB2 The pose of box 2 frame measured and expressed in frame F.
template <typename S>
void setUpBoxToBoxSimple(fcl::Box<S>* o1, fcl::Box<S>* o2,
                         fcl::Transform3<S>* X_FB1, fcl::Transform3<S>* X_FB2,
                         S delta_x = S(-0.6)) {
  const fcl::Vector3<S> box1_size(2, 2, 2);
  const fcl::Vector3<S> box2_size(1, 1, 2);
  // Box 1 is fixed.
  X_FB1->setIdentity();
  X_FB1->translation() << 0, 0, 1;
  X_FB2->setIdentity();
  X_FB2->translation() << delta_x, 0, 1;

  // Make the box
  fcl::Box<S> box1(box1_size);
  fcl::Box<S> box2(box2_size);
  *o1 = box1;
  *o2 = box2;
}

template <typename S>
void testBoxVsBoxSimple() {
  // Should collide for delta x less than 1.5
  const int test_n_intersect = 100;
  for (auto i = 0; i < test_n_intersect; i++) {
    S delta_x = 1.5 * double(i) / test_n_intersect;
    fcl::Box<S> b1, b2;
    fcl::Transform3<S> pose_1, pose_2;
    setUpBoxToBoxSimple(&b1, &b2, &pose_1, &pose_2, delta_x);

    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&b1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&b2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = pose_1.inverse() * pose_2;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto result = mpr.Intersect(minkowski_diff, &intersect_data);
    EXPECT_TRUE(mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data));
    EXPECT_TRUE(result == detail::MPR<S>::IntersectStatus::Intersect);
  }

  // Should separate for greater distance
  const int test_n_separate = 100;
  for (auto i = 0; i < test_n_separate; i++) {
    S delta_x = 1.5 + 1e-2 + 0.1 * (double(i) / test_n_separate);
    fcl::Box<S> b1, b2;
    fcl::Transform3<S> pose_1, pose_2;
    setUpBoxToBoxSimple(&b1, &b2, &pose_1, &pose_2, delta_x);

    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&b1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&b2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = pose_1.inverse() * pose_2;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    detail::MPR<S> mpr(1000, 1e-6);
    // auto result = mpr.Intersect(minkowski_diff);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto result = mpr.Intersect(minkowski_diff, &intersect_data);
    EXPECT_TRUE(result == detail::MPR<S>::IntersectStatus::Separated);
  }
}

template <typename S>
void testBoxVsBoxWithPrimitive() {
  fcl::Box<S> box1(0.1, 0.1, 0.1);
  fcl::Box<S> box2(0.2, 0.2, 0.2);
  int test_n = 1e6;
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> box1_pose, box2_pose;
    test::generateRandomTransform(extent, box1_pose);
    test::generateRandomTransform(extent, box2_pose);
    bool intersect_by_primitive =
        detail::boxBoxIntersect<S>(box1, box1_pose, box2, box2_pose, nullptr);

    // Make the minkowski difference
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&box2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = box1_pose.inverse() * box2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK/MPR
    detail::GJK<S> gjk(1000, 1e-6);
    auto gjk_result = gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (intersect_by_primitive) {
      if (gjk_result != detail::GJK<S>::Inside) {
        std::cout << "GJK failed to detect box collision " << std::endl;
      }
      EXPECT_TRUE(mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data));
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
    } else {
      if (gjk_result != detail::GJK<S>::Valid) {
        std::cout << "GJK failed to detect box separation " << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Separated);
    }
  }
}

template <typename S>
void testMPRWithCapsuleInstance() {
  fcl::Capsule<S> capsule_1(0.05, 0.2);
  fcl::Capsule<S> capsule_2(0.05, 0.2);
  fcl::Transform3<S> obj1_pose, obj2_pose;
  obj1_pose.setIdentity();
  obj2_pose.setIdentity();

  // Object 1 pose
  {
    Matrix3<S> rotation;
    rotation.setZero();
    rotation.row(0) = Vector3<S>(-0.5400293668847500062, -0.7492522363174671796,
                                 0.3833919264608088295);
    rotation.row(1) = Vector3<S>(0.06916196063564683527, 0.414480495836045626,
                                 0.9074263285647705679);
    rotation.row(2) = Vector3<S>(-0.838799681749515802, 0.5165530030352935009,
                                 -0.1720118860780580072);
    obj1_pose.linear().matrix() = rotation;
    obj1_pose.translation() = Vector3<S>(
        0.0453002929687499889, 0.2223999023437500111, 0.2935363769531250111);
  }

  // Object 2 pose
  {
    Matrix3<S> rotation;
    rotation.setZero();
    rotation.row(0) = Vector3<S>(-0.5389046582521762607, -0.7774713847560658087,
                                 0.3242221695066370146);
    rotation.row(1) = Vector3<S>(-0.01245611760117290578, 0.3922071179913819705,
                                 0.9197925971278524404);
    rotation.row(2) = Vector3<S>(-0.842274666868481181, 0.4916419657457848369,
                                 -0.2210465178852849544);
    obj2_pose.linear().matrix() = rotation;
    obj2_pose.translation() = Vector3<S>(
        0.044952392578125, 0.2228576660156250111, 0.2934448242187500111);
  }

  S distance;
  fcl::Vector3<S> point_1, point_2;
  detail::capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2, obj2_pose,
                                    &distance, &point_1, &point_2);
  const bool intersect_by_primitive = distance <= 0;
  std::cout << "Primitive intersect_by_primitive is " << intersect_by_primitive
            << std::endl;

  // Make the minkowski difference
  detail::MinkowskiDiff<S> minkowski_diff;
  minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
  minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
  minkowski_diff.support_function = detail::computeSupport<S>;
  minkowski_diff.interior_function = detail::computeInterior<S>;
  minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
  minkowski_diff.toshape1 =
      minkowski_diff.toshape0.inverse().rotation().matrix();

  // Try with GJK/MPR
  detail::GJK<S> gjk(1000, 1e-6);
  gjk.evaluate(minkowski_diff, Vector3<S>::UnitY());
  std::cout << gjk.getStringStatus() << std::endl;
  detail::MPR<S> mpr(1000, 1e-6);
  typename detail::MPR<S>::IntersectData intersect_data;
  auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
  if (mpr_result == detail::MPR<S>::IntersectStatus::Intersect) {
    std::cout << "MPR reports intersect" << std::endl;
    std::cout << minkowski_diff.interior().norm() << std::endl;
    const bool within_error_or_enclosed =
        std::abs(distance) < 1e-6 ||
        mpr.IsOriginEnclosed(minkowski_diff, intersect_data);
    EXPECT_TRUE(within_error_or_enclosed);
    std::cout << mpr.IsOriginEnclosed(minkowski_diff, intersect_data)
              << std::endl;
  } else if (mpr_result == detail::MPR<S>::IntersectStatus::Separated) {
    std::cout << "MPR reports seperated" << std::endl;
  } else {
    std::cout << "MPR reports failed" << std::endl;
  }
}

template <typename S>
void testMPRWithCapsulePair() {
  fcl::Capsule<S> capsule_1(0.05, 0.2);
  fcl::Capsule<S> capsule_2(0.05, 0.2);
  int test_n = 1e6;
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);
    S distance;
    fcl::Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2,
                                      obj2_pose, &distance, &point_1, &point_2);
    const bool intersect_by_primitive = distance <= 0;

    // Make the minkowski difference
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK/MPR
    detail::GJK<S> gjk(1000, 1e-6);
    auto gjk_result = gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (intersect_by_primitive) {
      if (gjk_result != detail::GJK<S>::Inside) {
        std::cout << "GJK failed to detect capsule collision" << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect ||
                  std::abs(distance) < 1e-6);
      const bool within_error_or_enclosed =
          std::abs(distance) < 1e-6 ||
          mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data);
      if (!within_error_or_enclosed) {
        std::cout << "Object1 pose" << std::endl;
        std::cout << obj1_pose.matrix() << std::endl;
        std::cout << "Object2 pose" << std::endl;
        std::cout << obj2_pose.matrix() << std::endl;
      }
      EXPECT_TRUE(within_error_or_enclosed);
    } else {
      if (gjk_result != detail::GJK<S>::Valid) {
        std::cout << "GJK failed to detect capsule separation " << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Separated);
    }
  }
}

// Callers
GTEST_TEST(MPR_Primitive_Test, box2box_simple_test) {
  testBoxVsBoxSimple<float>();
  testBoxVsBoxSimple<double>();
}

GTEST_TEST(MPR_Primitive_Test, box2box_with_primitive) {
  testBoxVsBoxWithPrimitive<float>();
  testBoxVsBoxWithPrimitive<double>();
}

GTEST_TEST(MPR_Primitive_Test, capsule2capsule_test) {
  testMPRWithCapsulePair<float>();
  testMPRWithCapsulePair<double>();
  // testMPRWithCapsuleInstance<double>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
