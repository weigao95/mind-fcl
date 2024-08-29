#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>

#include "fcl/cvx_collide/gjk.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_triangle.h"
#include "fcl/narrowphase/detail/shape_pair_intersect.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void testBoxWithTriangleP(const Box<S>& box, const TriangleP<S>& tri) {
  std::array<S, 6> xyz_extent{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  int test_n = 1e4;
  detail::GJKSolver<S> gjk_solver;
  detail::ShapePairIntersectSolver<S> intersect_solver(&gjk_solver);
  Simplex<S> triangle_simplex(tri.a, tri.b, tri.c);

  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> body_1_pose, body_2_pose;
    fcl::test::generateRandomTransform(xyz_extent, body_1_pose);
    fcl::test::generateRandomTransform(xyz_extent, body_2_pose);

    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&tri);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape1 =
        body_2_pose.linear().transpose() * body_1_pose.linear();
    minkowski_diff.toshape0 =
        body_1_pose.inverse(Eigen::Isometry) * body_2_pose;

    // Try with GJK/MPR
    detail::GJK2<S> gjk(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    cvx_collide::GJK_Status gjk_status = gjk.Evaluate(minkowski_diff, simplex);
    const bool is_overlap0 = (cvx_collide::GJK_Status::Intersect == gjk_status);

    // Run primitive with body_2_pose
    const bool is_overlap_test1 = detail::boxTriangleIntersect(
        box, body_1_pose, tri.a, tri.b, tri.c, body_2_pose);

    // Run solver without contact
    const bool is_overlap_test2 = gjk_solver.shapeTriangleIntersect(
        box, body_1_pose, tri.a, tri.b, tri.c, body_2_pose, nullptr);

    // Run solver with contact
    CollisionPenetrationContactData<S> penetration_contacts;
    const bool is_overlap_test3 =
        gjk_solver.shapeTriangleIntersect(box, body_1_pose, tri.a, tri.b, tri.c,
                                          body_2_pose, &penetration_contacts);

    // Apply with transformed point
    Vector3<S> transformed_a = body_2_pose * tri.a;
    Vector3<S> transformed_b = body_2_pose * tri.b;
    Vector3<S> transformed_c = body_2_pose * tri.c;

    // Run primitive without body_2_pose
    const bool is_overlap_test4 = detail::boxTriangleIntersect(
        box, body_1_pose, transformed_a, transformed_b, transformed_c);

    // Run solver without contact and tf2
    const bool is_overlap_test5 = gjk_solver.shapeTriangleIntersect(
        box, body_1_pose, transformed_a, transformed_b, transformed_c, nullptr);

    // Run solver without tf2, with contact
    penetration_contacts.clear();
    const bool is_overlap_test6 = gjk_solver.shapeTriangleIntersect(
        box, body_1_pose, transformed_a, transformed_b, transformed_c,
        &penetration_contacts);

    // Try intersect solver
    CollisionRequest<S> request;
    CollisionResult<S> result;
    // auto meta_fn = [](Contact<S>&) -> bool { return false; };
    detail::ContactMeta<S> contact_meta;
    intersect_solver.ShapeSimplexIntersect(box, body_1_pose, triangle_simplex,
                                           body_2_pose, request, contact_meta,
                                           result);
    const bool is_overlap_test7 = result.isCollision();

    request.useDefaultPenetration();
    result.clear();
    intersect_solver.ShapeSimplexIntersect(box, body_1_pose, triangle_simplex,
                                           body_2_pose, request, contact_meta,
                                           result);
    const bool is_overlap_test8 = result.isCollision();

    EXPECT_EQ(is_overlap_test1, is_overlap0);
    EXPECT_EQ(is_overlap_test2, is_overlap0);
    EXPECT_EQ(is_overlap_test3, is_overlap0);
    EXPECT_EQ(is_overlap_test4, is_overlap0);
    EXPECT_EQ(is_overlap_test5, is_overlap0);
    EXPECT_EQ(is_overlap_test6, is_overlap0);
    EXPECT_EQ(is_overlap_test7, is_overlap0);
    EXPECT_EQ(is_overlap_test8, is_overlap0);
  }
}

template <typename S>
void test_triangle_box() {
  const int test_n = 10;
  for (int i = 0; i < test_n; i++) {
    const auto a = test::generateRandomBox<S>();
    const auto b = test::generateRandomTriangleP<S>();
    testBoxWithTriangleP<S>(a, b);
  }
}

GTEST_TEST(FCL_Box_Triangle, collision_triangle_box) {
  test_triangle_box<float>();
  test_triangle_box<double>();
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
