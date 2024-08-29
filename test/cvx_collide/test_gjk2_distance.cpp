//
// Created by wei on 23-1-12.
//

#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/gjk.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename S>
void testSphereVsSphereSeparation() {
  const S sphere_radius = 1.0;
  const S expected_distance = 0.2;
  const S displacement = sphere_radius * 2 + expected_distance;

  fcl::Sphere<S> sphere1(sphere_radius);
  fcl::Sphere<S> sphere2(sphere_radius);
  int test_n = 1e6;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    test::generateRandomTransform(extent, sphere1_pose);
    test::generateRandomTransform(extent, sphere2_pose);
    sphere2_pose.translation() =
        sphere1_pose.translation() + direction_in_world * displacement;

    // Setup the env
    MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&sphere1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&sphere2);
    minkowski_diff.support_function = computeSupport<S>;
    minkowski_diff.interior_function = computeInterior<S>;
    minkowski_diff.toshape0 = sphere1_pose.inverse() * sphere2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run GJK2
    GJK2<S> gjk(1000, 1e-6);
    Vector3<S> guess = Vector3<S>::UnitX();
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_separation_output;
    GJK_Status gjk_result =
        gjk.Evaluate(minkowski_diff, simplex, guess, &min_separation_output);
    EXPECT_TRUE(gjk_result == GJK_Status::Separated);
    auto difference = std::abs(expected_distance -
                               min_separation_output.separation_distance());
    EXPECT_TRUE(difference < 1e-3);
    if (difference > 1e-3) {
      std::cout << "GJK2 failed to compute separation distance with error "
                << difference << std::endl;
      std::cout << "object1 pose" << std::endl;
      std::cout << sphere1_pose.matrix() << std::endl;
      std::cout << "object2 pose" << std::endl;
      std::cout << sphere2_pose.matrix() << std::endl;
    }
  }
}

template <typename S>
void testObjectPairInstance() {
  // Construct the object
  fcl::Capsule<S> object1(0.5, 1.3);
  fcl::Capsule<S> object2(0.5, 1.3);
  fcl::Transform3<S> object1_pose, object2_pose;

  // Object 1 pose
  {
    object1_pose.setIdentity();
    object1_pose.translation().x() = 1.8916015625;
    object1_pose.translation().y() = 1.1961669921875;
    object1_pose.translation().z() = 1.0020751953125;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        -0.1823781667231642312, 0.5834536857025155454, -0.7914038166087196124);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.4082009544663793843, 0.6873272332027842157, 0.6007938542216808564);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.8944887842088544705, 0.4326234750195090406, 0.112812870660265227);
    object1_pose.linear().matrix() = rotation_matrix;
  }

  // Object 2 pose
  {
    object2_pose.setIdentity();
    object2_pose.translation().x() = 1.6461181640625;
    object2_pose.translation().y() = 0.496826171875;
    object2_pose.translation().z() = 0.091796875;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.317527870519881128, -0.3295648710688394156, -0.8891361241117632375);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.1094271233639760865, 0.9186616351567586936, -0.3795872821412688003);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.9419138793517136676, 0.217825149715291877, 0.2556373369367547776);
    object2_pose.linear().matrix() = rotation_matrix;
  }

  // Setup the env
  MinkowskiDiff<S> minkowski_diff;
  minkowski_diff.shapes[0] = constructGJKGeometry(&object1);
  minkowski_diff.shapes[1] = constructGJKGeometry(&object2);
  minkowski_diff.support_function = computeSupport<S>;
  minkowski_diff.interior_function = computeInterior<S>;
  minkowski_diff.toshape0 = object1_pose.inverse() * object2_pose;
  minkowski_diff.toshape1 =
      minkowski_diff.toshape0.inverse().rotation().matrix();

  // Run GJK2
  GJK2<S> gjk(1000, 1e-7);
  Vector3<S> guess = Vector3<S>::UnitX();
  GJKSimplex<S> simplex;
  typename GJK2<S>::MinSeparationDistanceOutput min_separation_output;
  GJK_Status gjk_result =
      gjk.Evaluate(minkowski_diff, simplex, guess, &min_separation_output);
  EXPECT_TRUE(gjk_result == GJK_Status::Separated);
  std::cout << "min distance point norm is "
            << min_separation_output.separation_distance() << std::endl;
  std::cout << min_separation_output.point_on_minkowskidiff() << std::endl;

  // Check the error
  S capsule_distance;
  {
    Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(object1, object1_pose, object2,
                                      object2_pose, &capsule_distance, &point_1,
                                      &point_2);
  }

  auto difference =
      std::abs(min_separation_output.separation_distance() - capsule_distance);
  if (difference > 5e-3) {
    std::cout
        << "GJK2 failed to compute capsule separation distance with error "
        << min_separation_output.separation_distance() - capsule_distance
        << std::endl;
    std::cout << "Primitive algorithm produces " << capsule_distance
              << std::endl;
  }
}

template <typename S>
void testCapsuleVsCapsuleSeparation() {
  fcl::Capsule<S> capsule_1(0.5, 1.3);
  fcl::Capsule<S> capsule_2(0.5, 1.3);
  int test_n = 1e6;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};
  std::size_t separate_count = 0;
  std::size_t mismatch_count = 0;

  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);

    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&capsule_1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&capsule_2);
    minkowski_diff.support_function = computeSupport<S>;
    minkowski_diff.interior_function = computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run primitive
    S capsule_distance;
    fcl::Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2,
                                      obj2_pose, &capsule_distance, &point_1,
                                      &point_2);
    // Negative implies penetration
    if (capsule_distance <= 1e-6) continue;
    separate_count += 1;

    // Run GJK2
    GJK2<S> gjk(1000, 1e-6);
    Vector3<S> guess = Vector3<S>::UnitX();
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_separation_output;
    GJK_Status gjk_result =
        gjk.Evaluate(minkowski_diff, simplex, guess, &min_separation_output);
    EXPECT_TRUE(gjk_result == GJK_Status::Separated);

    // Check the error
    auto difference = std::abs(min_separation_output.separation_distance() -
                               capsule_distance);
    if (difference > 5e-3) {
      mismatch_count += 1;
      std::cout
          << "GJK2 failed to compute capsule separation distance with error "
          << min_separation_output.separation_distance() - capsule_distance
          << std::endl;
      std::cout << "object1 pose" << std::endl;
      std::cout << obj1_pose.matrix() << std::endl;
      std::cout << "object2 pose" << std::endl;
      std::cout << obj2_pose.matrix() << std::endl;
    }
  }

  // Logging
  EXPECT_TRUE(mismatch_count <= 2);
  std::cout << "Separate count is " << separate_count << std::endl;
  std::cout << "Mismatch count is " << mismatch_count << std::endl;
}

template <typename S>
void testConvexSphere() {
  const S sphere_radius = 1.0;
  const S expected_distance = 0.2;
  const S displacement = sphere_radius * 2 + expected_distance;

  auto sphere1 = std::make_shared<fcl::Convex<S>>(
      fcl::test::makeSphereConvex<S>(sphere_radius));
  auto sphere2 = std::make_shared<fcl::Convex<S>>(
      fcl::test::makeSphereConvex<S>(sphere_radius));
  int test_n = 1e5;
  int invalid_gjk = 0;
  int new_gjk_larger_distance_count = 0;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    test::generateRandomTransform(extent, sphere1_pose);
    test::generateRandomTransform(extent, sphere2_pose);
    sphere2_pose.translation() =
        sphere1_pose.translation() + direction_in_world * displacement;

    // Setup the env
    MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(sphere1.get());
    minkowski_diff.shapes[1] = constructGJKGeometry(sphere2.get());
    minkowski_diff.support_function = computeSupport<S>;
    minkowski_diff.interior_function = computeInterior<S>;
    minkowski_diff.toshape0 = sphere1_pose.inverse() * sphere2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run GJK2
    Vector3<S> guess = Vector3<S>::UnitX();
    S gjk2_min_distance = 0.0;
    GJK_Status gjk2_result;
    {
      GJK2<S> gjk(1000, 1e-6);
      GJKSimplex<S> simplex;
      typename GJK2<S>::MinSeparationDistanceOutput min_separation_output;
      gjk2_result =
          gjk.Evaluate(minkowski_diff, simplex, guess, &min_separation_output);
      gjk2_min_distance = min_separation_output.separation_distance();
    }
    if (gjk2_result != GJK_Status::Separated) {
      std::cout << "GJK2 failed to detect separation " << std::endl;
      std::cout << "object1 pose" << std::endl;
      std::cout << sphere1_pose.matrix() << std::endl;
      std::cout << "object2 pose" << std::endl;
      std::cout << sphere2_pose.matrix() << std::endl;
    }

    // Run GJK
    typename detail::GJK<S>::Status gjk_status;
    S gjk_min_distance = 0;
    {
      detail::GJK<S> gjk(1000, 1e-6);
      gjk_status = gjk.evaluate(minkowski_diff, guess);
      if (gjk_status == detail::GJK<S>::Valid) {
        Vector3<S> w0 = Vector3<S>::Zero();
        Vector3<S> w1 = Vector3<S>::Zero();
        for (size_t i = 0; i < gjk.getSimplex()->rank; ++i) {
          S p = gjk.getSimplex()->p[i];

          // Check the weight?
          if (p < -1e-2) {
            gjk_status = detail::GJK<S>::Failed;
          }

          // Use the weight to interpolate
          w0.noalias() +=
              minkowski_diff.support(gjk.getSimplex()->c[i]->d, 0) * p;
          w1.noalias() +=
              minkowski_diff.support(-gjk.getSimplex()->c[i]->d, 1) * p;
        }
        gjk_min_distance = (w0 - w1).norm();
      }
    }

    // Compare the result
    if (gjk_status == detail::GJK<S>::Valid) {
      S distance_diff = gjk2_min_distance - gjk_min_distance;
      if (distance_diff > 1e-3) {
        new_gjk_larger_distance_count += 1;
      }
    } else {
      invalid_gjk += 1;
    }
  }

  // Logging
  std::cout << "# of tests: " << test_n << std::endl;
  std::cout << "# of failed old gjk: " << invalid_gjk << std::endl;
  std::cout << "# of case where new gjk2 compute a larger distance: "
            << new_gjk_larger_distance_count << std::endl;
  EXPECT_TRUE(new_gjk_larger_distance_count <= 1000);
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(GJK2_Distance_Test, sphere_test) {
  fcl::detail::testSphereVsSphereSeparation<double>();
  fcl::detail::testSphereVsSphereSeparation<float>();

  // Debug
  // fcl::detail::testObjectPairInstance<double>();
}

GTEST_TEST(GJK2_Distance_Test, capsule_test) {
  fcl::detail::testCapsuleVsCapsuleSeparation<double>();

  // Debug method
  // fcl::detail::testObjectPairInstance<double>();
}

GTEST_TEST(GJK2_Distance_Test, convex_test) {
  fcl::detail::testConvexSphere<double>();

  // Debug
  // fcl::detail::testObjectPairInstance<double>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
