//
// Created by mech-mind_gw on 3/24/2022.
//

#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/gjk_shape.h"
#include "fcl/cvx_collide/mpr.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "retired_epa.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void testSphereVsSphere() {
  const S sphere_radius = 1.0;
  const S expected_penetration = 0.2;
  const S displacement = sphere_radius * 2 - expected_penetration;
  fcl::Sphere<S> sphere1(sphere_radius);
  fcl::Sphere<S> sphere2(sphere_radius);
  int test_n = 10000;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    test::generateRandomTransform(extent, sphere1_pose);
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    direction_in_world.normalize();
    sphere1_pose.setIdentity();
    sphere2_pose = sphere1_pose;
    sphere2_pose.translation() += direction_in_world * displacement;

    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&sphere1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&sphere2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = sphere1_pose.inverse() * sphere2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    Vector3<S> guess = Vector3<S>::UnitX();
    auto gjk_result = gjk.evaluate(minkowski_diff, guess);
    if (gjk_result == detail::GJK<S>::Inside) {
      detail::EPA<S> epa(128, 64, 255, 1e-6);
      typename detail::EPA<S>::Status epa_status = epa.evaluate(gjk, guess);
      if (epa_status != detail::EPA<S>::Failed) {
        if (std::abs(epa.depth - expected_penetration) > 1e-3) {
          // std::cout << "EPA depth is " << epa.depth << std::endl;
        }
      }

      // Run mpr
      detail::MPR<S> mpr(1000, 1e-6);
      typename detail::MPR<S>::DirectedPenetrationData data;
      fcl::Vector3<S> direction =
          sphere1_pose.inverse().linear().matrix() * direction_in_world;
      auto status = mpr.DirectedPenetration(minkowski_diff, direction, data);
      EXPECT_TRUE(status == detail::MPR<S>::DirectedPenetrationStatus::OK);
      EXPECT_TRUE(std::abs(data.distance_on_direction - expected_penetration) <
                  1e-5);
    } else {
      std::cout << "GJK failure" << std::endl;
    }
  }
}

template <typename S>
void testMoveToSeperated(const ShapeBase<S>* shape1, const ShapeBase<S>* shape2,
                         std::function<bool(const fcl::Transform3<S>&,
                                            const fcl::Transform3<S>&, double*)>
                             is_intersect_functor,
                         int test_n = 10,
                         double move_from_penetration_depth_eps = 1e-4,
                         double intersect_tolerance_upperbound = 1e-4) {
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  for (auto i = 0; i < test_n; i++) {
    fcl::Transform3<S> shape1_pose, shape2_pose;
    shape1_pose.setIdentity();
    test::generateRandomTransform(extent, shape2_pose);

    // Generate a random direction
    fcl::Vector3<S> direction_in_shape1;
    direction_in_shape1.setRandom();
    if (direction_in_shape1.norm() < 1e-6) {
      direction_in_shape1[0] += 0.1;
    }
    direction_in_shape1.normalize();

    // Make the diff
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(shape1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(shape2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = shape1_pose.inverse() * shape2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Determine the intersection
    if (!is_intersect_functor(shape1_pose, shape2_pose, nullptr)) {
      continue;
    }

    // If they intersect
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::DirectedPenetrationData data;
    const fcl::Vector3<S>& direction = direction_in_shape1;
    auto status = mpr.DirectedPenetration(minkowski_diff, direction, data);
    EXPECT_TRUE(status == detail::MPR<S>::DirectedPenetrationStatus::OK);

    // Move the shape2 by some distance to collision-free
    {
      fcl::Transform3<S> moved_shape2_pose = shape2_pose;
      moved_shape2_pose.translation() +=
          (data.distance_on_direction + move_from_penetration_depth_eps) *
          direction_in_shape1;

      // Check collision result
      double tolerance = 0.0;
      if (is_intersect_functor(shape1_pose, moved_shape2_pose, &tolerance)) {
        // It should be a small violation
        EXPECT_TRUE(tolerance < intersect_tolerance_upperbound);
        if (tolerance > intersect_tolerance_upperbound) {
          std::cout << "Wrong collision result after move with tolerance "
                    << tolerance << std::endl;
          std::cout << "The shape2 pose" << std::endl;
          std::cout << shape2_pose.matrix() << std::endl;
          std::cout << "The direction to test" << std::endl;
          std::cout << direction_in_shape1 << std::endl;
        }
      }
    }

    // Move the shape2 by some distance, but still collision
    {
      fcl::Transform3<S> moved_shape2_pose = shape2_pose;
      moved_shape2_pose.translation() +=
          (data.distance_on_direction - move_from_penetration_depth_eps) *
          direction_in_shape1;

      // Check collision result
      double tolerance = 0.0;
      if (!is_intersect_functor(shape1_pose, moved_shape2_pose, &tolerance)) {
        // It should be a small violation
        EXPECT_LE(tolerance, intersect_tolerance_upperbound);
        if (tolerance > intersect_tolerance_upperbound) {
          std::cout << "Wrong collision-free result after move with tolerance "
                    << tolerance << " " << intersect_tolerance_upperbound << std::endl;
        }
      }
    }
  }
}

template <typename S>
void testShapeVsShapeInstance(const ShapeBase<S>* shape1,
                              const ShapeBase<S>* shape2) {
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};

  fcl::Transform3<S> shape1_pose, shape2_pose;
  shape1_pose.setIdentity();
  shape2_pose.setIdentity();
  // shape2_pose.translation()[1] += S(1.8);
  shape2_pose.translation().x() = 0.0472595214843750111;
  shape2_pose.translation().y() = -0.09856567382812500555;
  shape2_pose.translation().z() = 0.0498229980468750111;

  // The rotation matrix
  fcl::Matrix3<S> rotation_matrix;
  rotation_matrix.row(0) = fcl::Vector3<S>(
      0.001102700129736739441, 0.001327474267536466216, -0.9999985109311377851);
  rotation_matrix.row(1) = fcl::Vector3<S>(
      0.9835928429985736354, -0.1804006789176861303, 0.0008451319549329973108);
  rotation_matrix.row(2) = fcl::Vector3<S>(
      -0.1803992883977296735, -0.9835923102882142555, -0.001504622040860770444);
  shape2_pose.linear().matrix() = rotation_matrix;

  // The direction to test
  Vector3<S> direction = Vector3<S>::UnitY();
  direction.normalize();

  // Setup the env
  detail::MinkowskiDiff<S> minkowski_diff;
  minkowski_diff.shapes[0] = constructGJKGeometry(shape1);
  minkowski_diff.shapes[1] = constructGJKGeometry(shape2);
  minkowski_diff.support_function = detail::computeSupport<S>;
  minkowski_diff.interior_function = detail::computeInterior<S>;
  minkowski_diff.toshape0 = shape1_pose.inverse() * shape2_pose;
  minkowski_diff.toshape1 =
      minkowski_diff.toshape0.inverse().rotation().matrix();

  // Run GJK/EPA
  detail::GJK<S> gjk(1000, 1e-6);
  Vector3<S> guess = Vector3<S>::UnitY();
  auto gjk_result = gjk.evaluate(minkowski_diff, guess);
  if (gjk_result == detail::GJK<S>::Inside) {
    detail::EPA<S> epa(128, 64, 255, 1e-6);
    typename detail::EPA<S>::Status epa_status = epa.evaluate(gjk, guess);
    if (epa_status != detail::EPA<S>::Failed) {
      std::cout << "The epa depth " << epa.depth << std::endl;

      detail::MPR<S> mpr(1000, 1e-6);
      typename detail::MPR<S>::DirectedPenetrationData data;
      auto status = mpr.DirectedPenetration(minkowski_diff, direction, data);
      EXPECT_TRUE(status == detail::MPR<S>::DirectedPenetrationStatus::OK);
      std::cout << "MPR depth " << data.distance_on_direction << std::endl;
    }
  }
}

template <typename S>
void sphereMoveToSeparateTest() {
  fcl::Sphere<S> sphere1(S(0.15));
  fcl::Sphere<S> sphere2(S(0.15));
  auto intersect_functor = [](const fcl::Transform3<S>& pose1,
                              const fcl::Transform3<S>& pose2,
                              double* tolerance) -> bool {
    fcl::Vector3<S> center_diff = pose1.translation() - pose2.translation();
    if (tolerance != nullptr) *tolerance = center_diff.norm() - 0.3;
    return center_diff.norm() < 0.3;
  };

  // Start testing
  int test_n = 1e5;
  fcl::testMoveToSeperated<S>(&sphere1, &sphere2, intersect_functor, test_n);
}

template <typename S>
void boxMoveToSeparateTest() {
  fcl::Box<S> box1(0.2, 0.2, 0.2);
  fcl::Box<S> box2(0.2, 0.2, 0.2);
  auto intersect_functor = [&box1, &box2](const fcl::Transform3<S>& pose1,
                                          const fcl::Transform3<S>& pose2,
                                          double* tolerance) -> bool {
    fcl::CollisionPenetrationContactData<S> contacts;
    int return_code;
    fcl::Vector3<S> normal;
    S depth;
    bool intersect_by_primitive =
        fcl::detail::boxBox2<S>(box1.side, pose1, box2.side, pose2, normal,
                                &depth, &return_code, 4, contacts);
    if (tolerance != nullptr) {
      *tolerance = depth;
    }
    return intersect_by_primitive;
  };

  int test_n = 1e5;
  fcl::testMoveToSeperated<S>(&box1, &box2, intersect_functor, test_n);
}

template <typename S>
void convexSphereMoveToSeparateTest() {
  S radius = 0.15;
  auto sphere1 = fcl::test::makeSphereConvex<S>(radius);
  auto sphere2 = fcl::test::makeSphereConvex<S>(radius);
  auto intersect_functor = [&sphere1, &sphere2, &radius](
                               const fcl::Transform3<S>& pose1,
                               const fcl::Transform3<S>& pose2,
                               double* tolerance) -> bool {
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&sphere1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&sphere2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = pose1.inverse() * pose2;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    detail::MPR<S> mpr(1000, 1e-8);
    auto mpr_result = mpr.Intersect(minkowski_diff, nullptr);
    bool is_intersect =
        (mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
    if (tolerance != nullptr) {
      fcl::Vector3<S> center_diff = pose1.translation() - pose2.translation();
      *tolerance = center_diff.norm() - (2.0 * radius);
    }
    return is_intersect;
  };

  // Start testing
  int test_n = 1e5;
  fcl::testMoveToSeperated<S>(&sphere1, &sphere2, intersect_functor, test_n);
}

}  // namespace fcl

GTEST_TEST(MPR_Penetration_Test, sphere_test) {
  fcl::testSphereVsSphere<double>();
  fcl::testSphereVsSphere<float>();
}

GTEST_TEST(MPR_Penetration_Test, sphere_move_to_separate) {
  fcl::sphereMoveToSeparateTest<double>();
  fcl::sphereMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_Penetration_Test, box_move_to_separate) {
  fcl::boxMoveToSeparateTest<double>();
  fcl::boxMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_Penetration_Test, sphere_convex_move_to_separate) {
  fcl::convexSphereMoveToSeparateTest<double>();
  fcl::convexSphereMoveToSeparateTest<float>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
