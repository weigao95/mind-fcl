//
// Created by mech-mind_gw on 2/10/2022.
//

#include <gtest/gtest.h>

#include <memory>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/mpr.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename T>
void sphereConvexTest(bool use_nv_sphere = false) {
  T radius = 0.1;
  std::shared_ptr<fcl::Convex<T>> sphere_1, sphere_2;
  if (use_nv_sphere) {
    sphere_1 = std::make_shared<fcl::Convex<T>>(
        fcl::test::makeSphereConvex<T>(radius));
    sphere_2 = std::make_shared<fcl::Convex<T>>(
        fcl::test::makeSphereConvex<T>(radius));
  } else {
    sphere_1 = std::make_shared<fcl::Convex<T>>(
        fcl::test::makeSphereConvexUV<T>(radius, 20, 20));
    sphere_2 = std::make_shared<fcl::Convex<T>>(
        fcl::test::makeSphereConvexUV<T>(radius, 20, 20));
  }

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
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(sphere_1.get());
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(sphere_2.get());
    minkowski_diff.support_function = detail::computeSupport<T>;
    minkowski_diff.interior_function = detail::computeInterior<T>;
    minkowski_diff.toshape0 = body_1_pose.inverse() * body_2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK and MPR
    detail::GJK<T> gjk(1000, 1e-6);
    auto gjk_result = gjk.evaluate(minkowski_diff, Vector3<T>::UnitX());
    detail::MPR<T> mpr(1000, 1e-6);
    typename detail::MPR<T>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);

    if (gjk_result == detail::GJK<T>::Inside) {
      if (mpr_result != detail::MPR<T>::IntersectStatus::Intersect) {
        std::cout << "*******************************" << std::endl;
        std::cout << "MPR mismatch with GJK, which claims intersect"
                  << std::endl;
        auto distance =
            (body_1_pose.translation() - body_2_pose.translation()).norm();
        std::cout << "Sphere primitive is_collide " << (distance <= 2 * radius)
                  << " with Gt distance " << distance << " radius " << radius
                  << std::endl;
        std::cout << "MPR reports "
                  << fcl::test::intersectStatusToString<T>(mpr_result)
                  << std::endl;
        std::cout << "MPR PortalEncloseOrigin "
                  << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
                  << std::endl;
      }
    }

    if (gjk_result == detail::GJK<T>::Valid) {
      if (mpr_result != detail::MPR<T>::IntersectStatus::Separated) {
        std::cout << "*******************************" << std::endl;
        std::cout << "MPR mismatch with GJK, which claims seperated with "
                     "min-distance to be "
                  << gjk.distance << std::endl;
        auto distance =
            (body_1_pose.translation() - body_2_pose.translation()).norm();
        std::cout << "Sphere primitive is_collide " << (distance <= 2 * radius)
                  << " with Gt distance " << distance << " radius " << radius
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
GTEST_TEST(MPR_Convex_Test, test_sphere_convex) {
  sphereConvexTest<float>(false);
  sphereConvexTest<double>(false);
}

GTEST_TEST(MPR_Convex_Test, test_sphere_convex_uv) {
  sphereConvexTest<float>(true);
  sphereConvexTest<double>(true);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
