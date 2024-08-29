//
// Created by mech-mind_gw on 6/8/2022.
//

#include <gtest/gtest.h>

#include <memory>

#include "fcl/cvx_collide/epa.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "test_fcl_utility.h"
#include "create_primitive_mesh.h"
#include "retired_epa.h"
#include "retired_gjk.h"

namespace fcl {
namespace detail {

template <typename S>
void testSphereVsSpherePenetration() {
  const S sphere_radius = 1.0;
  const S expected_penetration = 0.2;
  const S displacement = sphere_radius * 2 - expected_penetration;
  fcl::Sphere<S> sphere1(sphere_radius);
  fcl::Sphere<S> sphere2(sphere_radius);
  int test_n = 100000;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    direction_in_world.normalize();
    sphere1_pose.setIdentity();
    sphere2_pose = sphere1_pose;
    sphere2_pose.translation() += direction_in_world * displacement;

    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&sphere1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&sphere2);
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
      EPA2<S> epa;
      S depth;
      // auto epa_status =
      //     epa.Evaluate(gjk, minkowski_diff, &depth, nullptr, nullptr);
      auto epa_status = runEPA2_WithOldGJK(epa, gjk, minkowski_diff, &depth);
      if (epa_status != EPA_Status::Failed) {
        if (std::abs(depth - expected_penetration) > 1e-4) {
          std::cout << "EPA depth wrong with difference "
                    << depth - expected_penetration << std::endl;
        }
      } else {
        std::cout << "EPA failure" << std::endl;
      }
    } else {
      std::cout << "GJK failure" << std::endl;
    }
  }
}

template <typename S>
void testBoxVsBoxPenetration() {
  fcl::Box<S> box1(2, 1, 0.5);
  fcl::Box<S> box2(1, 1, 1);
  int test_n = 1000000;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};
  std::size_t collide_count = 0;

  // Info about mismatch
  std::size_t mismatch_count = 0;
  S max_mismatch = 0.0;
  fcl::Transform3<S> obj1_pose_mismatch, obj2_pose_mismatch;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);

    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&box1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&box2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run primitive
    CollisionPenetrationContactData<S> contacts;
    int return_code;
    fcl::Vector3<S> normal;
    S depth;
    bool intersect_by_primitive =
        fcl::detail::boxBox2<S>(box1.side, obj1_pose, box2.side, obj2_pose,
                                normal, &depth, &return_code, 4, contacts);
    if (!intersect_by_primitive) continue;
    collide_count += 1;

    // NOTE(wei): debug
    // std::cout << "Running " << test_idx << std::endl;
    // std::cout << obj1_pose.matrix() << std::endl;
    // std::cout << obj2_pose.matrix() << std::endl;

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    EPA2<S> epa;
    Vector3<S> guess = Vector3<S>::UnitX();
    auto gjk_result = gjk.evaluate(minkowski_diff, guess);
    if (gjk_result == detail::GJK<S>::Inside) {
      S depth_epa;
      // auto epa_status =
      //     epa.Evaluate(gjk, minkowski_diff, &depth_epa, nullptr, nullptr);
      auto epa_status = runEPA2_WithOldGJK(epa, gjk, minkowski_diff, &depth_epa);
      if (epa_status != EPA_Status::Failed) {
        const S mismatch_depth = std::abs(depth - depth_epa);
        if (mismatch_depth > 1e-4 && depth < depth_epa) {
          mismatch_count += 1;
          if (mismatch_depth > max_mismatch) {
            max_mismatch = mismatch_depth;
            obj1_pose_mismatch = obj1_pose;
            obj2_pose_mismatch = obj2_pose;
          }
        }
      } else {
        std::cout << "EPA failure" << std::endl;
        std::cout << obj1_pose.matrix() << std::endl;
        std::cout << obj2_pose.matrix() << std::endl;
      }
    } else {
      // std::cout << "GJK failure" << std::endl;
    }
  }

  // Logging
  std::cout << "Collision count is " << collide_count << std::endl;
  if (mismatch_count > 0) {
    std::cout << "Mismatch count is " << mismatch_count << " with max "
              << max_mismatch << std::endl;
    std::cout << obj1_pose_mismatch.matrix() << std::endl;
    std::cout << obj2_pose_mismatch.matrix() << std::endl;
  }
}

template <typename S>
void testCapsuleVsCapsulePenetration() {
  fcl::Capsule<S> capsule_1(0.5, 1.3);
  fcl::Capsule<S> capsule_2(0.5, 1.3);
  int test_n = 1000000;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};
  std::size_t collide_count = 0;

  // Info about mismatch
  std::size_t mismatch_count = 0;
  S max_mismatch = 0.0;
  S max_mismatch_gt_depth = 0.0;
  fcl::Transform3<S> obj1_pose_mismatch, obj2_pose_mismatch;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);

    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&capsule_1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&capsule_2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
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
    if (capsule_distance >= 0) continue;
    capsule_distance *= -1;
    collide_count += 1;

    // NOTE(wei): debug
    // std::cout << "Running " << test_idx << std::endl;
    // std::cout << obj1_pose.matrix() << std::endl;
    // std::cout << obj2_pose.matrix() << std::endl;

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    Vector3<S> guess = Vector3<S>::UnitX();
    auto gjk_result = gjk.evaluate(minkowski_diff, guess);
    if (gjk_result == detail::GJK<S>::Inside) {
      EPA2<S> epa(512);
      S depth_epa;
      // auto epa_status =
      //     epa.Evaluate(gjk, minkowski_diff, &depth_epa, nullptr, nullptr);
      auto epa_status = runEPA2_WithOldGJK(epa, gjk, minkowski_diff, &depth_epa);
      if (epa_status != EPA_Status::Failed) {
        const S mismatch_depth = std::abs(capsule_distance - depth_epa);
        if (mismatch_depth > 1e-4 && capsule_distance < depth_epa) {
          mismatch_count += 1;
          if (mismatch_depth > max_mismatch) {
            max_mismatch = mismatch_depth;
            obj1_pose_mismatch = obj1_pose;
            obj2_pose_mismatch = obj2_pose;
            max_mismatch_gt_depth = capsule_distance;
          }
        }
      } else {
        // std::cout << "EPA failure" << std::endl;
      }
    } else {
      std::cout << "GJK failure" << std::endl;
    }
  }

  // Logging
  std::cout << "Collision count is " << collide_count << std::endl;
  if (mismatch_count > 0) {
    std::cout << "Mismatch count is " << mismatch_count << " with max "
              << max_mismatch << std::endl;
    std::cout << obj1_pose_mismatch.matrix() << std::endl;
    std::cout << obj2_pose_mismatch.matrix() << std::endl;
    std::cout << "The GT is " << max_mismatch_gt_depth << std::endl;
  }
}

template <typename S>
void testPenetrationInstance() {
  // fcl::Box<S> object1(2, 1, 0.5);
  // fcl::Box<S> object2(1, 1, 1);
  fcl::Capsule<S> object1(0.5, 1.3);
  fcl::Capsule<S> object2(0.5, 1.3);

  // The object pose
  fcl::Transform3<S> obj1_pose, obj2_pose;

  // obj1 pose
  {
    auto& shape1_pose = obj1_pose;
    shape1_pose.setIdentity();
    shape1_pose.translation().x() = -0.01806640625;
    shape1_pose.translation().y() = -1.669677734375;
    shape1_pose.translation().z() = -1.6229248046875;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.507635951042175293, 0.7651183009147644043, 0.3961057364940643311);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        0.6602797508239746094, -0.6408261656761169434, 0.391627967357635498);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.5534766316413879395, 0.0627361685037612915, -0.8304985165596008301);
    shape1_pose.linear().matrix() = rotation_matrix;
  }

  // obj2 pose
  {
    auto& shape2_pose = obj2_pose;
    shape2_pose.setIdentity();
    shape2_pose.translation().x() = 0.6041259765625;
    shape2_pose.translation().y() = -0.8148193359375;
    shape2_pose.translation().z() = -1.3558349609375;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.1389396339654922485, 0.990229487419128418, 0.0118880709633231163);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        0.8491457104682922363, -0.125303804874420166, 0.5130793452262878418);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.5095559358596801758, -0.06119235977530479431, -0.8582587838172912598);
    shape2_pose.linear().matrix() = rotation_matrix;
  }

  // The test code
  {
    // Setup the env
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = constructGJKGeometry(&object1);
    minkowski_diff.shapes[1] = constructGJKGeometry(&object2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    Vector3<S> guess = Vector3<S>::UnitX();
    auto gjk_result = gjk.evaluate(minkowski_diff, guess);
    if (gjk_result == detail::GJK<S>::Inside) {
      EPA2<S> epa(512);
      S depth_epa;
      Vector3<S> p1, p2;
      // auto epa_status = epa.Evaluate(gjk, minkowski_diff, &depth_epa, &p1, &p2);
      auto epa_status = runEPA2_WithOldGJK(epa, gjk, minkowski_diff, &depth_epa, &p1, &p2);

      // EPA works
      if (epa_status != EPA_Status::Failed) {
        std::cout << "EPA depth is " << depth_epa << std::endl;
      } else {
        std::cout << "EPA failure" << std::endl;
      }
    } else {
      std::cout << "GJK failure" << std::endl;
    }
  }

  {
    // Run primitive
    S capsule_distance;
    fcl::Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(object1, obj1_pose, object2, obj2_pose,
                                      &capsule_distance, &point_1, &point_2);
    std::cout << "Actual signed distance " << -capsule_distance << std::endl;
  }
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(EPA2_Test, sphere_test) {
  fcl::detail::testSphereVsSpherePenetration<double>();
  fcl::detail::testSphereVsSpherePenetration<float>();
}

GTEST_TEST(EPA2_Test, box_test) {
  fcl::detail::testBoxVsBoxPenetration<double>();
  fcl::detail::testBoxVsBoxPenetration<float>();
  // fcl::detail::testPenetrationInstance<double>();
}

GTEST_TEST(EPA2_Test, capsule_test) {
  fcl::detail::testCapsuleVsCapsulePenetration<double>();
  fcl::detail::testCapsuleVsCapsulePenetration<float>();
  // fcl::detail::testPenetrationInstance<float>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
