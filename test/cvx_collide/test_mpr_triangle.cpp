//
// Created by wei on 2/8/22.
//

#include <gtest/gtest.h>

#include <memory>

#include "fcl/cvx_collide/mpr.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void testBoxVsTriangleNoTransform() {
  fcl::Box<S> box(0.1, 0.1, 0.1);
  int test_n = 1e5;
  for (auto test_i = 0; test_i < test_n; test_i++) {
    // A random triangle
    Vector3<S> a, b, c;
    a.setRandom();
    b.setRandom();
    c.setRandom();
    fcl::TriangleP<S> triangle(a, b, c);

    // Make the minkowski difference
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangle);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0.setIdentity();
    minkowski_diff.toshape1.setIdentity();

    // Try with GJK and MPR
    detail::GJK<S> gjk(1000, 1e-6);
    auto gjk_result = gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (gjk_result == detail::GJK<S>::Inside) {
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
      if (mpr_result != detail::MPR<S>::IntersectStatus::Intersect) {
        std::cout << "A triangle should intersect " << triangle << std::endl;
      }
    }

    if (gjk_result == detail::GJK<S>::Valid) {
      // Should be no intersection
      if (mpr_result != detail::MPR<S>::IntersectStatus::Separated) {
        std::cout << "A triangle that is seperated from GJK " << triangle
                  << " with min-distance to be " << gjk.distance << std::endl;
        std::cout << "MPR reports "
                  << fcl::test::intersectStatusToString<S>(mpr_result)
                  << std::endl;
        std::cout << "MPR PortalEncloseOrigin "
                  << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
                  << std::endl;
      }
    }
  }
}

template <typename S>
void testBoxVsTriangleWithTransform() {
  fcl::Box<S> box(0.1, 0.1, 0.1);
  int test_n = 1e5;
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  for (auto i = 0; i < test_n; i++) {
    fcl::Transform3<S> box_pose, triangle_pose;
    test::generateRandomTransform(extent, box_pose);
    test::generateRandomTransform(extent, triangle_pose);

    // A random triangle
    Vector3<S> a, b, c;
    a.setRandom();
    b.setRandom();
    c.setRandom();
    fcl::TriangleP<S> triangle(a, b, c);

    // Make the minkowski difference
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangle);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = box_pose.inverse() * triangle_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();

    // Try with GJK and MPR
    detail::GJK<S> gjk(1000, 1e-6);
    auto gjk_result = gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    detail::MPR<S> mpr(1000, 1e-6);
    typename detail::MPR<S>::IntersectData intersect_data;
    auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
    if (gjk_result == detail::GJK<S>::Inside) {
      EXPECT_TRUE(mpr_result == detail::MPR<S>::IntersectStatus::Intersect);
      if (mpr_result != detail::MPR<S>::IntersectStatus::Intersect) {
        std::cout << "A triangle should intersect " << triangle << std::endl;
      }
    }

    if (gjk_result == detail::GJK<S>::Valid) {
      // Should be no intersection
      if (mpr_result != detail::MPR<S>::IntersectStatus::Separated) {
        std::cout << "A triangle that is seperated from GJK " << triangle
                  << " with min-distance to be " << gjk.distance << std::endl;
        std::cout << "MPR reports "
                  << fcl::test::intersectStatusToString<S>(mpr_result)
                  << std::endl;
        std::cout << "MPR PortalEncloseOrigin "
                  << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
                  << std::endl;
      }
    }
  }
}

template <typename S>
void triangleBoxInstanceTest() {
  fcl::Box<S> box(0.1, 0.1, 0.1);
  fcl::Vector3<S> a(-0.452864, -0.57152, -0.149449);
  fcl::Vector3<S> b(-0.3961, 0.631947, -0.356487);
  fcl::Vector3<S> c(0.442061, -0.81402, 0.616321);
  fcl::TriangleP<S> triangle(a, b, c);

  // Make the minkowski difference
  detail::MinkowskiDiff<S> minkowski_diff;
  minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
  minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangle);
  minkowski_diff.support_function = detail::computeSupport<S>;
  minkowski_diff.interior_function = detail::computeInterior<S>;
  minkowski_diff.toshape0.setIdentity();
  minkowski_diff.toshape1.setIdentity();

  // Try with GJK and MPR
  detail::GJK<S> gjk(1000, 1e-6);
  gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
  detail::MPR<S> mpr(1000, 1e-6);
  typename detail::MPR<S>::IntersectData intersect_data;
  auto mpr_result = mpr.Intersect(minkowski_diff, &intersect_data);
  std::cout << "GJK status is " << gjk.getStringStatus() << " with distance "
            << gjk.distance << std::endl;
  std::cout << "MPR reports "
            << fcl::test::intersectStatusToString<S>(mpr_result) << std::endl;
  std::cout << "MPR PortalEncloseOrigin "
            << mpr.IsOriginEnclosedDebug(minkowski_diff, intersect_data)
            << std::endl;
}

// Callers
GTEST_TEST(MPR_Test, box2triangle_no_transform_test) {
  triangleBoxInstanceTest<float>();
  triangleBoxInstanceTest<double>();
  testBoxVsTriangleNoTransform<double>();
}

GTEST_TEST(MPR_Test, box2triangle_transfomed_test) {
  testBoxVsTriangleWithTransform<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
