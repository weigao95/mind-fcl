//
// Created by mech-mind_gw on 3/28/2022.
//

#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/gjk_shape.h"
#include "fcl/cvx_collide/mpr.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_capsule.h"
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
  int test_n = 100000;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    // test::generateRandomTransform(extent, sphere1_pose);
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    direction_in_world.normalize();
    // fcl::Vector3<S> direction_in_world = Vector3<S>::UnitX();
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

      // Change the direction a bit
      fcl::Vector3<S> direction =
          sphere1_pose.inverse().linear().matrix() * direction_in_world;
      fcl::Vector3<S> direction_disturb;
      direction_disturb.setRandom();
      direction_disturb *= 0.3;
      fcl::Vector3<S> disturbed_direction = direction + direction_disturb;
      if (disturbed_direction.squaredNorm() <= 0) {
        disturbed_direction[0] += 0.3;
      }
      disturbed_direction.normalize();

      // Run mpr
      detail::MPR<S> mpr(1000, 1e-6);
      typename detail::MPR<S>::IncrementalMinimumPenetrationData data;
      auto status = mpr.IncrementalMinimumPenetrationDistance(
          minkowski_diff, disturbed_direction, data);
      // Debug
      if (status != detail::MPR<S>::IncrementalPenetrationStatus::OK) {
        std::cout << "Distributed direction" << std::endl;
        std::cout << disturbed_direction << std::endl;
        std::cout << "tf1" << std::endl;
        std::cout << sphere1_pose.matrix() << std::endl;
        std::cout << "tf2" << std::endl;
        std::cout << sphere2_pose.matrix() << std::endl;
        std::cout << "****************" << std::endl;
      }
      // End of debug
      EXPECT_TRUE(status == detail::MPR<S>::IncrementalPenetrationStatus::OK);
      EXPECT_LE(std::abs(data.minimum_penetration - expected_penetration),
                S(2e-5));
    } else {
      std::cout << "GJK failure" << std::endl;
    }
  }
}

template <typename S>
void testRefinementMoveToSeperated(
    const ShapeBase<S>* shape1, const ShapeBase<S>* shape2,
    std::function<bool(const fcl::Transform3<S>&, const fcl::Transform3<S>&,
                       fcl::Vector3<S>* direction, S* depth)>
        is_intersect_functor,
    int test_n = 10, S disturb_angle_tan = 0.577,
    S move_from_penetration_depth_eps = 1e-4,
    S intersect_tolerance_upperbound = 1e-5) {
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  for (auto test_i = 0; test_i < test_n; test_i++) {
    fcl::Transform3<S> shape1_pose, shape2_pose;
    shape1_pose.setIdentity();
    test::generateRandomTransform(extent, shape2_pose);

    // Determine the intersection
    S gt_depth;
    Vector3<S> direction;
    bool is_intersect =
        is_intersect_functor(shape1_pose, shape2_pose, &direction, &gt_depth);
    if (!is_intersect) {
      continue;
    }

    // Make the shape
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(shape1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(shape2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = shape1_pose.inverse() * shape2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    Vector3<S> guess = direction;
    auto gjk_result = gjk.evaluate(minkowski_diff, guess);
    if (gjk_result == detail::GJK<S>::Inside) {
      detail::EPA<S> epa(256, 256, 1000, 1e-6);
      typename detail::EPA<S>::Status epa_status = epa.evaluate(gjk, guess);
      if (epa_status != detail::EPA<S>::Failed) {
        // std::cout << "The epa depth " << epa.depth << std::endl;
      }
    }

    // Change the direction a bit
    Vector3<S> disturbed_direction;
    test::disturbDirectionByAngle(direction, disturb_angle_tan,
                                  disturbed_direction);

    // Run mpr
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IncrementalMinimumPenetrationData data;
    auto status = mpr.IncrementalMinimumPenetrationDistance(
        minkowski_diff, disturbed_direction, data);
    EXPECT_TRUE(status == detail::MPR<S>::IncrementalPenetrationStatus::OK);
    if (status != detail::MPR<S>::IncrementalPenetrationStatus::OK) {
      std::cout << "The shape2 pose" << std::endl;
      std::cout << shape2_pose.matrix() << std::endl;
      std::cout << "The disturbed_direction direction" << std::endl;
      std::cout << disturbed_direction << std::endl;
    }

    // Test the distance
    if (data.minimum_penetration > gt_depth + move_from_penetration_depth_eps) {
      std::cout << "Possible non-optimal minimum penetration computation"
                << std::endl;
      std::cout << "Ground-truth depth " << gt_depth << std::endl;
      std::cout << "MPR depth " << data.minimum_penetration << std::endl;
      /*std::cout << "The shape2 pose" << std::endl;
      std::cout << shape2_pose.matrix() << std::endl;
      std::cout << "The disturbed direction" << std::endl;
      std::cout << disturbed_direction << std::endl;
      std::cout << "The initial direction" << std::endl;
      std::cout << direction << std::endl;*/
    }

    // Check that the direction/distance lead to separation
    {
      fcl::Transform3<S> moved_shape2_pose = shape2_pose;
      moved_shape2_pose.translation() +=
          (data.minimum_penetration + move_from_penetration_depth_eps) *
          data.penetration_direction;

      // Check collision result
      S tolerance = 0.0;
      if (is_intersect_functor(shape1_pose, moved_shape2_pose, nullptr,
                               &tolerance)) {
        // It should be a small violation
        EXPECT_TRUE(tolerance < intersect_tolerance_upperbound);
        /*if (tolerance > intersect_tolerance_upperbound) {
          std::cout << "Wrong collision result after move with tolerance " <<
        tolerance << std::endl; std::cout << "The shape2 pose" << std::endl;
          std::cout << shape2_pose.matrix() << std::endl;
          std::cout << "The computed direction" << std::endl;
          std::cout << data.penetration_direction << std::endl;
        }*/
      }
    }

    // Check that moving not sufficient distance lead to penetration
    {
      fcl::Transform3<S> moved_shape2_pose = shape2_pose;
      moved_shape2_pose.translation() +=
          (data.minimum_penetration - move_from_penetration_depth_eps) *
          data.penetration_direction;

      // Check collision result
      S tolerance = 0.0;
      if (!is_intersect_functor(shape1_pose, moved_shape2_pose, nullptr,
                                &tolerance)) {
        // It should be a small violation
        EXPECT_TRUE(tolerance < intersect_tolerance_upperbound);
        if (tolerance > intersect_tolerance_upperbound) {
          std::cout << "Wrong collision-free result after move with tolerance "
                    << tolerance << std::endl;
        }
      }
    }
  }
}

template <typename S>
void testShapeVsShapeInstance(const ShapeBase<S>* shape1,
                              const ShapeBase<S>* shape2) {
  fcl::Transform3<S> shape1_pose, shape2_pose;
  shape1_pose.setIdentity();
  shape2_pose.setIdentity();
  shape2_pose.translation().x() = 0.6013973355293273926;
  shape2_pose.translation().y() = -1.496152520179748535;
  shape2_pose.translation().z() = 0.7999053001403808594;

  // The rotation matrix
  fcl::Matrix3<S> rotation_matrix;
  rotation_matrix.row(0) = fcl::Vector3<S>(
      -0.3656753808078196943, -0.2969471611190888094, -0.8821019778769175756);
  rotation_matrix.row(1) = fcl::Vector3<S>(
      -0.5573794451611232548, -0.6891388675721128454, 0.4630505105421053313);
  rotation_matrix.row(2) = fcl::Vector3<S>(
      -0.7453922926774428914, 0.6609916827803054007, 0.08648887392217365078);
  // shape2_pose.linear().matrix() = rotation_matrix;

  // The direction to test
  Vector3<S> direction = Vector3<S>::UnitY();
  direction.normalize();

  // Setup the env
  detail::MinkowskiDiff<S> minkowski_diff;
  minkowski_diff.shapes[0] = detail::constructGJKGeometry(shape1);
  minkowski_diff.shapes[1] = detail::constructGJKGeometry(shape2);
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
    }
  }

  // Run mpr
  Vector3<S> mpr_direction_in(0.5068867802619934082, -0.8381253480911254883,
                              0.2015235275030136108);
  detail::MPR<S> mpr(1000, 1e-6);
  typename detail::MPR<S>::IncrementalMinimumPenetrationData data;
  auto status = mpr.IncrementalMinimumPenetrationDistance(
      minkowski_diff, mpr_direction_in, data);
  EXPECT_TRUE(status == detail::MPR<S>::IncrementalPenetrationStatus::OK);
  std::cout << "MPR depth " << data.minimum_penetration << std::endl;
}

template <typename S>
void sphereRefinementMoveToSeparateTest() {
  fcl::Sphere<S> sphere1(S(0.15));
  fcl::Sphere<S> sphere2(S(0.15));
  auto intersect_functor = [](const fcl::Transform3<S>& pose1,
                              const fcl::Transform3<S>& pose2,
                              Vector3<S>* direction, S* depth) -> bool {
    fcl::Vector3<S> center_diff = pose2.translation() - pose1.translation();
    if (direction != nullptr) {
      *direction = center_diff.normalized();
    }
    if (depth != nullptr) {
      *depth = 0.3 - center_diff.norm();
    }
    return center_diff.norm() < 0.3;
  };

  // Start testing
  // fcl::testShapeVsShapeInstance(&sphere1, &sphere2);
  int test_n = 10000;
  fcl::testRefinementMoveToSeperated<S>(&sphere1, &sphere2, intersect_functor,
                                        test_n);
}

template <typename S>
void sphereConvexRefinementMoveToSeparateTest() {
  S radius = 0.15;
  auto sphere1 = fcl::test::makeSphereConvex<S>(radius);
  auto sphere2 = fcl::test::makeSphereConvex<S>(radius);
  auto intersect_functor = [&sphere1, &sphere2](const fcl::Transform3<S>& pose1,
                                                const fcl::Transform3<S>& pose2,
                                                Vector3<S>* direction,
                                                S* depth) -> bool {
    // Rough estimation on distance
    fcl::Vector3<S> center_diff = pose2.translation() - pose1.translation();
    if (direction != nullptr) {
      *direction = center_diff.normalized();
    }
    if (depth != nullptr) {
      *depth = 0.3 - center_diff.norm();
    }

    // Make the shape
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&sphere1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&sphere2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = pose1.inverse() * pose2;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Intersect is by mpr
    detail::MPR<S> mpr(1000, 1e-8);
    auto mpr_result = mpr.Intersect(minkowski_diff, nullptr);
    bool is_intersect =
        (mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
    return is_intersect;
  };

  // Start testing
  int test_n = 10000;
  S disturb_angle_tan = 0.2;
  fcl::testRefinementMoveToSeperated<S>(&sphere1, &sphere2, intersect_functor,
                                        test_n, disturb_angle_tan);
}

template <typename S>
void capsuleRefinementMoveToSeparateTest() {
  fcl::Capsule<S> capsule_1(0.15, 0.3);
  fcl::Capsule<S> capsule_2(0.15, 0.3);
  auto intersect_functor = [&capsule_1, &capsule_2](
                               const fcl::Transform3<S>& pose1,
                               const fcl::Transform3<S>& pose2,
                               Vector3<S>* direction, S* depth) -> bool {
    S capsule_distance;
    fcl::Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(capsule_1, pose1, capsule_2, pose2,
                                      &capsule_distance, &point_1, &point_2);
    if (depth != nullptr) {
      *depth = -capsule_distance;
    }
    if (direction != nullptr) {
      *direction = point_1 - point_2;
      direction->normalize();
    }
    return capsule_distance <= 0;
  };

  // Start testing
  int test_n = 10000;
  fcl::testRefinementMoveToSeperated<S>(&capsule_1, &capsule_2,
                                        intersect_functor, test_n);
}

template <typename S>
void sphereCapsuleRefinementMoveToSeparateTest() {
  fcl::Sphere<S> sphere_1(0.2);
  fcl::Capsule<S> capsule_2(0.15, 0.3);
  auto intersect_functor = [&sphere_1, &capsule_2](
                               const fcl::Transform3<S>& pose1,
                               const fcl::Transform3<S>& pose2,
                               Vector3<S>* direction, S* depth) -> bool {
    fcl::Vector3<S> point_1, point_2;
    CollisionPenetrationContactData<S> contact;
    bool intersect = detail::sphereCapsuleIntersect<S>(
        sphere_1, pose1, capsule_2, pose2, &contact);
    if (intersect && depth != nullptr) {
      assert(contact.size() == 1);
      *depth = contact[0].penetration_depth;
    }
    if (intersect && direction != nullptr) {
      *direction = contact[0].normal;
      direction->normalize();
    }
    return intersect;
  };

  // Start testing
  int test_n = 10000;
  fcl::testRefinementMoveToSeperated<S>(&sphere_1, &capsule_2,
                                        intersect_functor, test_n);
}

template <typename S>
void boxRefinementMoveToSeparateTest() {
  fcl::Box<S> box1(0.2, 0.2, 0.2);
  fcl::Box<S> box2(0.2, 0.2, 0.2);
  auto intersect_functor = [&box1, &box2](const fcl::Transform3<S>& pose1,
                                          const fcl::Transform3<S>& pose2,
                                          Vector3<S>* direction,
                                          S* penetration) -> bool {
    fcl::CollisionPenetrationContactData<S> contacts;
    int return_code;
    fcl::Vector3<S> normal;
    S depth;
    bool intersect_by_primitive =
        fcl::detail::boxBox2<S>(box1.side, pose1, box2.side, pose2, normal,
                                &depth, &return_code, 4, contacts);
    if (penetration != nullptr) {
      *penetration = depth;
    }
    if (direction != nullptr) {
      *direction = normal;
    }
    return intersect_by_primitive;
  };

  // Start testing
  // fcl::testShapeVsShapeInstance(&box1, &box2);
  int test_n = 10000;
  S disturb_angle_tan = 0.2;
  fcl::testRefinementMoveToSeperated<S>(&box1, &box2, intersect_functor, test_n,
                                        disturb_angle_tan);
}

}  // namespace fcl

GTEST_TEST(MPR_LocalRefine_Test, sphere_test) {
  fcl::testSphereVsSphere<double>();
  fcl::testSphereVsSphere<float>();
}

GTEST_TEST(MPR_LocalRefine_Test, sphere_move_test) {
  fcl::sphereRefinementMoveToSeparateTest<double>();
  fcl::sphereRefinementMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_LocalRefine_Test, sphere_convex_move_test) {
  fcl::sphereConvexRefinementMoveToSeparateTest<double>();
  fcl::sphereConvexRefinementMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_LocalRefine_Test, capsule_move_test) {
  fcl::capsuleRefinementMoveToSeparateTest<double>();
  fcl::capsuleRefinementMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_LocalRefine_Test, sphere_capsule_move_test) {
  fcl::sphereCapsuleRefinementMoveToSeparateTest<double>();
  fcl::sphereCapsuleRefinementMoveToSeparateTest<float>();
}

GTEST_TEST(MPR_LocalRefine_Test, box_move_test) {
  fcl::boxRefinementMoveToSeparateTest<double>();
  fcl::boxRefinementMoveToSeparateTest<float>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
