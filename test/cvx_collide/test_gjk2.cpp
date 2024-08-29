//
// Created by wei on 22-6-18.
//
#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/epa.h"
#include "fcl/cvx_collide/gjk.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename S>
void testGJK2BoxVsBoxWithPrimitive() {
  fcl::Box<S> box1(0.1, 0.1, 0.1);
  fcl::Box<S> box2(0.2, 0.2, 0.2);
  int test_n = 1e6;
  std::array<S, 6> extent{-0.2, -0.2, -0.2, 0.2, 0.2, 0.2};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> box1_pose, box2_pose;
    test::generateRandomTransform(extent, box1_pose);
    test::generateRandomTransform(extent, box2_pose);
    int return_code = 0;
    CollisionPenetrationContactData<S> contacts;
    Vector3<S> normal;
    S depth = 0;
    detail::boxBox2(box1.side, box1_pose, box2.side, box2_pose, normal, &depth,
                    &return_code, 4, contacts);
    bool intersect_by_primitive = (return_code != 0);

    // Make the minkowski difference
    MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&box1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&box2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = box1_pose.inverse() * box2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK/MPR
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    auto gjk_status = gjk.Evaluate(minkowski_diff, simplex);
    MPR<S> mpr(1000, 1e-6);
    typename MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (intersect_by_primitive) {
      if (gjk_status != GJK_Status::Intersect) {
        std::cout << "GJK failed to detect box collision at depth " << depth
                  << std::endl;
        std::cout << "Box1 pose" << std::endl;
        std::cout << box1_pose.matrix() << std::endl;
        std::cout << "Box2 pose" << std::endl;
        std::cout << box2_pose.matrix() << std::endl;
      }
      EXPECT_TRUE(gjk_status == GJK_Status::Intersect);
      // EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
    } else {
      if (gjk_status != GJK_Status::Separated) {
        std::cout << "GJK failed to detect box separation at depth " << depth
                  << std::endl;
        std::cout << "Box1 pose" << std::endl;
        std::cout << box1_pose.matrix() << std::endl;
        std::cout << "Box2 pose" << std::endl;
        std::cout << box2_pose.matrix() << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Separated);
    }
  }
}

template <typename S>
void testGJK2WithCapsulePair() {
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
    capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2, obj2_pose,
                              &distance, &point_1, &point_2);
    const bool intersect_by_primitive = distance <= 0;

    // Make the minkowski difference
    MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&capsule_1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&capsule_2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK/MPR
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    auto gjk_result = gjk.Evaluate(minkowski_diff, simplex);
    MPR<S> mpr(1000, 1e-6);
    typename MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (intersect_by_primitive) {
      if (gjk_result != GJK_Status::Intersect) {
        std::cout << "GJK failed to detect capsule collision at distance "
                  << distance << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect ||
                  std::abs(distance) < 1e-6);
    } else {
      if (gjk_result != GJK_Status::Separated) {
        std::cout << "GJK failed to detect capsule separation at distance "
                  << distance << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Separated);
    }
  }
}

template <typename S>
void testGJK2BoxVsBoxWithPrimitiveInstance() {
  fcl::Box<S> obj1(0.1, 0.1, 0.1);
  fcl::Box<S> obj2(0.2, 0.2, 0.2);
  fcl::Transform3<S> obj1_pose, obj2_pose;
  obj1_pose.setIdentity();
  obj2_pose.setIdentity();

  // obj1 pose
  {
    auto& shape1_pose = obj1_pose;
    shape1_pose.setIdentity();
    shape1_pose.translation().x() = 0.1331561010330915562;
    shape1_pose.translation().y() = -0.04835435003042221069;
    shape1_pose.translation().z() = -0.08343856278806925653;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        -0.6001922922342994848, -0.6452300220595174052, 0.4727022646186623822);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.6144479604359862623, 0.750289910702751861, 0.2439646568946272354);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        -0.5120770608595877071, -0.1440252357426682894, -0.8467784924114895029);
    shape1_pose.linear().matrix() = rotation_matrix;
  }

  // obj2 pose
  {
    auto& shape2_pose = obj2_pose;
    shape2_pose.setIdentity();
    shape2_pose.translation().x() = -0.00932837892323731821;
    shape2_pose.translation().y() = 0.08853187374770643547;
    shape2_pose.translation().z() = -0.05020954832434654236;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.7155820814515156947, 0.2218768144779907359, -0.6623541076362519098);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        0.5764451852026616363, -0.7231199678065238778, 0.3805370686492457466);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        -0.3945291284077582228, -0.6541163438996084878, -0.6453515131160874052);
    shape2_pose.linear().matrix() = rotation_matrix;
  }

  // Try pose
  {
    /*bool intersect_by_primitive = detail::boxBoxIntersect<S>(
        obj1, obj1_pose, obj2, obj2_pose, nullptr);*/
    int return_code;
    CollisionPenetrationContactData<S> contacts;
    Vector3<S> normal;
    S depth;
    detail::boxBox2(obj1.side, obj1_pose, obj2.side, obj2_pose, normal, &depth,
                    &return_code, 4, contacts);
    bool intersect_by_primitive = (return_code != 0);

    // Make the minkowski difference
    MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&obj1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&obj2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK/MPR
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    auto gjk_status = gjk.Evaluate(minkowski_diff, simplex);
    MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (intersect_by_primitive) {
      if (gjk_status != GJK_Status::Intersect) {
        std::cout << "GJK failed to detect box collision at depth " << depth
                  << std::endl;
      }
      EXPECT_TRUE(gjk_status == GJK_Status::Intersect);
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
    } else {
      if (gjk_status != GJK_Status::Separated) {
        std::cout << "GJK failed to detect box separation at depth " << depth
                  << std::endl;
        std::cout << "Box1 pose" << std::endl;
        std::cout << obj1_pose.matrix() << std::endl;
        std::cout << "Box2 pose" << std::endl;
        std::cout << obj2_pose.matrix() << std::endl;
      }
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Separated);
    }
  }
}

template <typename T>
void sphereConvexTestGJK2ComparedMPR() {
  T radius = 0.1;
  auto sphere_1 =
      std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));
  auto sphere_2 =
      std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));

  // The extent used to generate the pose
  T extent = radius * 2;
  std::array<T, 6> xyz_extent{-extent, -extent, -extent,
                              extent,  extent,  extent};
  int test_n = 1e6;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<T> body_1_pose, body_2_pose;
    fcl::test::generateRandomTransform(xyz_extent, body_1_pose);
    fcl::test::generateRandomTransform(xyz_extent, body_2_pose);

    // Make the minkowski difference
    detail::MinkowskiDiff<T> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(sphere_1.get());
    minkowski_diff.shapes[1] = constructGJKGeometry(sphere_2.get());
    minkowski_diff.support_function = detail::computeSupport<T>;
    minkowski_diff.interior_function = detail::computeInterior<T>;
    minkowski_diff.toshape0 = body_1_pose.inverse() * body_2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK and MPR
    detail::GJK2<T> gjk(1000, 1e-6);
    GJKSimplex<T> simplex;
    auto gjk_result = gjk.Evaluate(minkowski_diff, simplex);
    detail::MPR<T> mpr(1000, 1e-6);
    typename detail::MPR<T>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);

    if (gjk_result == GJK_Status::Intersect) {
      if (mpr_result != detail::MPR<T>::IntersectStatus::Intersect) {
        std::cout << "MPR mismatch with GJK, which claims intersect"
                  << std::endl;
        std::cout << "MPR reports "
                  << fcl::test::intersectStatusToString<T>(mpr_result)
                  << std::endl;
        std::cout << "MPR PortalEncloseOrigin "
                  << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
                  << std::endl;
      }
    }

    if (gjk_result == GJK_Status::Separated) {
      if (mpr_result != detail::MPR<T>::IntersectStatus::Separated) {
        std::cout << "MPR mismatch with GJK, which claims seperated "
                  << std::endl;
        std::cout << "MPR reports "
                  << fcl::test::intersectStatusToString<T>(mpr_result)
                  << std::endl;
        std::cout << "MPR PortalEncloseOrigin "
                  << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
                  << std::endl;
      }
    }
  }
}

// Callers
GTEST_TEST(GJK2_Primitive_Test, box2box_with_primitive) {
  testGJK2BoxVsBoxWithPrimitiveInstance<double>();
  testGJK2BoxVsBoxWithPrimitive<double>();
  testGJK2BoxVsBoxWithPrimitive<float>();
}

GTEST_TEST(GJK2_Primitive_Test, capsule2capsule_test) {
  testGJK2WithCapsulePair<double>();
  testGJK2WithCapsulePair<float>();
}

GTEST_TEST(MPR_Convex_Test, test_sphere_convex_with_mpr) {
  sphereConvexTestGJK2ComparedMPR<float>();
  sphereConvexTestGJK2ComparedMPR<double>();
}

}  // namespace detail
}  // namespace fcl

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
