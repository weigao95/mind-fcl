//
// Created by mech-mind_gw on 1/10/2023.
//

#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/epa.h"
#include "fcl/cvx_collide/gjk.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "retired_epa.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace detail {

template <typename S>
void testSpherePenetration() {
  const S sphere_radius = 1.0;
  const S expected_penetration = 0.2;
  const S displacement = sphere_radius * 2 - expected_penetration;
  fcl::Sphere<S> sphere1(sphere_radius);
  fcl::Sphere<S> sphere2(sphere_radius);
  int test_n = 1e4;
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

    // Run GJK2
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    EPA2<S> epa;
    auto gjk_status = gjk.Evaluate(minkowski_diff, simplex);

    // Run EPA2
    if (gjk_status == GJK_Status::Intersect) {
      S depth;
      auto epa_status =
          epa.Evaluate(simplex, minkowski_diff, &depth, nullptr, nullptr);
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
void testBoxPenetration(bool use_cached_polytope) {
  fcl::Box<S> box1(2, 1, 0.5);
  fcl::Box<S> box2(1, 1, 1);
  int test_n = 1e6;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};
  std::size_t collide_count = 0;

  // Info about mismatch
  int mismatch_count = 0;
  S max_mismatch = 0.0;
  EPA2<S> epa;
  cvx_collide::Polytope<S> cached_polytope(epa.max_n_faces());
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

    // Run GJK2
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    auto gjk_status = gjk.Evaluate(minkowski_diff, simplex);
    if (gjk_status == GJK_Status::Intersect) {
      S depth_epa;
      EPA_Status epa_status;
      if (!use_cached_polytope) {
        epa_status =
            epa.Evaluate(simplex, minkowski_diff, &depth_epa, nullptr, nullptr);
      } else {
        epa_status = epa.Evaluate(simplex, minkowski_diff, &depth_epa, nullptr,
                                  nullptr, &cached_polytope);
      }
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
      std::cout << "GJK failure" << std::endl;
      std::cout << obj1_pose.matrix() << std::endl;
      std::cout << obj2_pose.matrix() << std::endl;
    }
  }

  // Logging
  std::cout << "Collision count is " << collide_count << std::endl;
  EXPECT_EQ(mismatch_count, 0);
  if (mismatch_count > 0) {
    std::cout << "Mismatch count is " << mismatch_count << " with max "
              << max_mismatch << std::endl;
    std::cout << obj1_pose_mismatch.matrix() << std::endl;
    std::cout << obj2_pose_mismatch.matrix() << std::endl;
    // std::cout << "The GT is " << max_mismatch_gt_depth << std::endl;
  }
}

template <typename S>
void testCapsulePenetration(bool use_cached_polytope) {
  fcl::Capsule<S> capsule_1(0.5, 1.3);
  fcl::Capsule<S> capsule_2(0.5, 1.3);
  int test_n = 1e6;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};
  std::size_t collide_count = 0;
  std::size_t max_n_faces = 512;

  // Info about mismatch
  std::size_t mismatch_count = 0;
  S max_mismatch = 0.0;
  S max_mismatch_gt_depth = 0.0;
  fcl::Transform3<S> obj1_pose_mismatch, obj2_pose_mismatch;
  cvx_collide::Polytope<S> cached_polytope(max_n_faces);
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
    Vector3<S> guess = Vector3<S>::UnitX();
    GJK2<S> gjk(1000, 1e-6);
    GJKSimplex<S> simplex;
    EPA2<S> epa(512);
    auto gjk_status = gjk.Evaluate(minkowski_diff, simplex, guess);
    if (gjk_status == GJK_Status::Intersect) {
      S depth_epa;
      EPA_Status epa_status;
      if (!use_cached_polytope) {
        epa_status =
            epa.Evaluate(simplex, minkowski_diff, &depth_epa, nullptr, nullptr);
      } else {
        epa_status = epa.Evaluate(simplex, minkowski_diff, &depth_epa, nullptr,
                                  nullptr, &cached_polytope);
      }

      // Check the result
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
        std::cout << "EPA failure" << std::endl;
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

}  // namespace detail
}  // namespace fcl

GTEST_TEST(EPA2_W_GJK2_Test, sphere_test) {
  fcl::detail::testSpherePenetration<double>();
  fcl::detail::testSpherePenetration<float>();
}

GTEST_TEST(EPA2_W_GJK2_Test, box_test) {
  fcl::detail::testBoxPenetration<double>(false);
  fcl::detail::testBoxPenetration<float>(false);
}

GTEST_TEST(EPA2_W_GJK2_Test, box_cached_polytope_test) {
  fcl::detail::testBoxPenetration<double>(true);
  fcl::detail::testBoxPenetration<float>(true);
}

GTEST_TEST(EPA2_W_GJK2_Test, capsule_test) {
  fcl::detail::testCapsulePenetration<double>(false);
  fcl::detail::testCapsulePenetration<float>(false);
}

GTEST_TEST(EPA2_W_GJK2_Test, capsule_cached_polytope_test) {
  fcl::detail::testCapsulePenetration<double>(true);
  fcl::detail::testCapsulePenetration<float>(true);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
