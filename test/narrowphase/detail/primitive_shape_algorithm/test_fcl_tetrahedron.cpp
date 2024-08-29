#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>

#include "fcl/geometry/shape/tetrahedron.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/tetrahedron_intersect.h"
#include "fcl/narrowphase/detail/shape_pair_intersect.h"
#include "test_fcl_utility.h"

using namespace fcl;

/// @brief Function to convert Tetrahedron to Convex shape
template <typename S>
fcl::Convex<S> tetrahedronToConvex(const Tetrahedron<S>& tetrahedron) {
  const std::vector<Vector3<S>> vertices(tetrahedron.vertices.begin(),
                                         tetrahedron.vertices.end());
  std::vector<int> faces;
  faces.reserve(16);
  for (int i = 0; i < 4; i++) {
    faces.emplace_back(3);
    faces.emplace_back((i + 1) % 4);
    faces.emplace_back((i + 2) % 4);
    faces.emplace_back((i + 3) % 4);
  }

  fcl::Convex<S> convex(
      std::make_shared<const std::vector<Vector3<S>>>(vertices), 4,
      std::make_shared<const std::vector<int>>(std::move(faces)));
  convex.computeLocalAABB();
  return convex;
}

template <typename S>
void testTetrahedronWithBox(const Tetrahedron<S>& tetrahedron,
                            const Box<S>& box) {
  const auto tetrahedron_convex = tetrahedronToConvex<S>(tetrahedron);
  Simplex<S> tetrahedron_simplex(
      tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2],
      tetrahedron.vertices[3]);
  std::array<S, 6> xyz_extent{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  int test_n = 1e4;
  detail::GJKSolver<S> gjk_solver;
  detail::ShapePairIntersectSolver<S> intersect_solver(&gjk_solver);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> tetrahedron_pose, box_pose;
    fcl::test::generateRandomTransform(xyz_extent, tetrahedron_pose);
    fcl::test::generateRandomTransform(xyz_extent, box_pose);
    const bool is_overlap = gjk_solver.shapeIntersect(
        box, box_pose, tetrahedron_convex, tetrahedron_pose, nullptr);

    const bool is_overlap_test1 = gjk_solver.shapeIntersect(
        box, box_pose, tetrahedron, tetrahedron_pose, nullptr);
    const bool is_overlap_test2 = detail::boxTerahedronIntersect(
        box, box_pose, tetrahedron, tetrahedron_pose);
    CollisionPenetrationContactData<S> penetration_contacts;
    const bool is_overlap_test3 = gjk_solver.shapeIntersect(
        box, box_pose, tetrahedron, tetrahedron_pose, &penetration_contacts);

    // Try intersect solver
    CollisionRequest<S> request;
    CollisionResult<S> result;
    // auto meta_fn = [](Contact<S>&) -> bool { return false; };
    detail::ContactMeta<S> contact_meta;
    intersect_solver.ShapeSimplexIntersect(box, box_pose, tetrahedron_simplex,
                                           tetrahedron_pose, request,
                                           contact_meta, result);
    const bool is_overlap_test4 = result.isCollision();

    request.useDefaultPenetration();
    result.clear();
    intersect_solver.ShapeSimplexIntersect(box, box_pose, tetrahedron_simplex,
                                           tetrahedron_pose, request,
                                           contact_meta, result);
    const bool is_overlap_test5 = result.isCollision();

    EXPECT_EQ(is_overlap_test1, is_overlap);
    EXPECT_EQ(is_overlap_test2, is_overlap);
    EXPECT_EQ(is_overlap_test3, is_overlap);
    EXPECT_EQ(is_overlap_test4, is_overlap);
    EXPECT_EQ(is_overlap_test5, is_overlap);
  }
}

template <typename S>
void testTetrahedronWithTriangleP(const Tetrahedron<S>& tetrahedron,
                                  const TriangleP<S>& triangle) {
  const fcl::Convex<S> tetrahedron_convex = tetrahedronToConvex<S>(tetrahedron);
  Simplex<S> tetrahedron_simplex(
      tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2],
      tetrahedron.vertices[3]);
  Simplex<S> triangle_simplex(triangle.a, triangle.b, triangle.c);
  detail::GJKSolver<S> gjk_solver;
  detail::ShapePairIntersectSolver<S> intersect_solver(&gjk_solver);

  std::array<S, 6> xyz_extent{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  int test_n = 1e4;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> triangle_pose;
    fcl::Transform3<S> tetrahedron_pose;
    fcl::test::generateRandomTransform(xyz_extent, triangle_pose);
    fcl::test::generateRandomTransform(xyz_extent, tetrahedron_pose);

    // Shape triangle interface with convex
    const bool is_overlap = gjk_solver.shapeTriangleIntersect(
        tetrahedron_convex, tetrahedron_pose, triangle.a, triangle.b,
        triangle.c, triangle_pose);

    // Primitive
    const bool is_overlap_test1 = detail::triangleTerahedronIntersect(
        triangle, triangle_pose, tetrahedron, tetrahedron_pose);

    // Using the shape interface
    const bool is_overlap_test2 = gjk_solver.shapeIntersect(
        triangle, triangle_pose, tetrahedron, tetrahedron_pose, nullptr);
    CollisionPenetrationContactData<S> penetration_contacts;
    const bool is_overlap_test3 =
        gjk_solver.shapeIntersect(triangle, triangle_pose, tetrahedron,
                                  tetrahedron_pose, &penetration_contacts);

    // Shape interface by reversed
    const bool is_overlap_test4 = gjk_solver.shapeIntersect(
        tetrahedron, tetrahedron_pose, triangle, triangle_pose, nullptr);
    penetration_contacts.clear();
    const bool is_overlap_test5 =
        gjk_solver.shapeIntersect(tetrahedron, tetrahedron_pose, triangle,
                                  triangle_pose, &penetration_contacts);

    EXPECT_EQ(is_overlap_test1, is_overlap);
    EXPECT_EQ(is_overlap_test2, is_overlap);
    EXPECT_EQ(is_overlap_test3, is_overlap);
    EXPECT_EQ(is_overlap_test4, is_overlap);
    EXPECT_EQ(is_overlap_test5, is_overlap);

    // Try with other solver
    {
      const auto& s1 = triangle_simplex;
      const auto& tf1 = triangle_pose;
      const auto& s2 = tetrahedron_simplex;
      const auto& tf2 = tetrahedron_pose;

      // Compute the relative transform
      Matrix3<S> rotation_2to1;
      Vector3<S> translation_2in1;
      relativeTransform(tf1.linear(), tf1.translation(), tf2.linear(),
                        tf2.translation(), rotation_2to1, translation_2in1);

      // Try intersect solver
      CollisionRequest<S> request;
      CollisionResult<S> result;
      // auto meta_fn = [](Contact<S>&) -> bool { return false; };
      detail::ContactMeta<S> contact_meta;
      intersect_solver.SimplexIntersect(s1, tf1, s2, tf2, rotation_2to1,
                                        translation_2in1, request, contact_meta,
                                        result);
      const bool is_overlap_test7 = result.isCollision();

      request.useDefaultPenetration();
      result.clear();
      intersect_solver.SimplexIntersect(s1, tf1, s2, tf2, rotation_2to1,
                                        translation_2in1, request, contact_meta,
                                        result);
      const bool is_overlap_test8 = result.isCollision();
      EXPECT_EQ(is_overlap_test7, is_overlap);
      EXPECT_EQ(is_overlap_test8, is_overlap);
    }

    // Reverse and try again
    {
      const auto& s1 = tetrahedron_simplex;
      const auto& tf1 = tetrahedron_pose;
      const auto& s2 = triangle_simplex;
      const auto& tf2 = triangle_pose;

      // Compute the relative transform
      Matrix3<S> rotation_2to1;
      Vector3<S> translation_2in1;
      relativeTransform(tf1.linear(), tf1.translation(), tf2.linear(),
                        tf2.translation(), rotation_2to1, translation_2in1);

      // Try intersect solver
      CollisionRequest<S> request;
      CollisionResult<S> result;
      // auto meta_fn = [](Contact<S>&) -> bool { return false; };
      detail::ContactMeta<S> contact_meta;
      intersect_solver.SimplexIntersect(s1, tf1, s2, tf2, rotation_2to1,
                                        translation_2in1, request, contact_meta,
                                        result);
      const bool is_overlap_test7 = result.isCollision();

      request.useDefaultPenetration();
      result.clear();
      intersect_solver.SimplexIntersect(s1, tf1, s2, tf2, rotation_2to1,
                                        translation_2in1, request, contact_meta,
                                        result);
      const bool is_overlap_test8 = result.isCollision();
      EXPECT_EQ(is_overlap_test7, is_overlap);
      EXPECT_EQ(is_overlap_test8, is_overlap);
    }
  }
}

template <typename S>
void testTetrahedronWithTetrahedron(const Tetrahedron<S>& a,
                                    const Tetrahedron<S>& b) {
  const auto a_convex = tetrahedronToConvex<S>(a);
  const auto b_convex = tetrahedronToConvex<S>(b);
  Simplex<S> a_simplex(a.vertices[0], a.vertices[1], a.vertices[2],
                       a.vertices[3]);
  Simplex<S> b_simplex(b.vertices[0], b.vertices[1], b.vertices[2],
                       b.vertices[3]);

  std::array<S, 6> xyz_extent{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  int test_n = 1e4;
  detail::GJKSolver<S> gjk_solver;
  detail::ShapePairIntersectSolver<S> intersect_solver(&gjk_solver);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> a_pose, b_pose;
    fcl::test::generateRandomTransform(xyz_extent, a_pose);
    fcl::test::generateRandomTransform(xyz_extent, b_pose);

    // GT by convex
    const bool is_overlap =
        gjk_solver.shapeIntersect(a_convex, a_pose, b_convex, b_pose, nullptr);

    // Primitive
    const bool is_overlap_test1 =
        detail::tetrahedronTrahedronIntersect(a, a_pose, b, b_pose);

    // GJK solve with/without contact
    const bool is_overlap_test2 =
        gjk_solver.shapeIntersect(a, a_pose, b, b_pose, nullptr);
    CollisionPenetrationContactData<S> penetration_contacts;
    const bool is_overlap_test3 =
        gjk_solver.shapeIntersect(a, a_pose, b, b_pose, &penetration_contacts);

    // GJK solver with tetrahedron
    const bool is_overlap_test4 = gjk_solver.shapeTetrahedronIntersect(
        a, a_pose, b.vertices[0], b.vertices[1], b.vertices[2], b.vertices[3],
        b_pose, nullptr);
    penetration_contacts.clear();
    const bool is_overlap_test5 = gjk_solver.shapeTetrahedronIntersect(
        a, a_pose, b.vertices[0], b.vertices[1], b.vertices[2], b.vertices[3],
        b_pose, &penetration_contacts);

    // Reverse of above
    const bool is_overlap_test6 = gjk_solver.shapeTetrahedronIntersect(
        b, b_pose, a.vertices[0], a.vertices[1], a.vertices[2], a.vertices[3],
        a_pose, nullptr);
    penetration_contacts.clear();
    const bool is_overlap_test7 = gjk_solver.shapeTetrahedronIntersect(
        b, b_pose, a.vertices[0], a.vertices[1], a.vertices[2], a.vertices[3],
        a_pose, &penetration_contacts);

    // Compute the relative transform
    Matrix3<S> rotation_b_to_a;
    Vector3<S> translation_b_in_a;
    relativeTransform(a_pose.linear(), a_pose.translation(), b_pose.linear(),
                      b_pose.translation(), rotation_b_to_a,
                      translation_b_in_a);

    // Try intersect solver
    CollisionRequest<S> request;
    CollisionResult<S> result;
    // auto meta_fn = [](Contact<S>&) -> bool { return false; };
    detail::ContactMeta<S> contact_meta;
    intersect_solver.SimplexIntersect(a_simplex, a_pose, b_simplex, b_pose,
                                      rotation_b_to_a, translation_b_in_a,
                                      request, contact_meta, result);
    const bool is_overlap_test8 = result.isCollision();

    request.useDefaultPenetration();
    result.clear();
    intersect_solver.SimplexIntersect(a_simplex, a_pose, b_simplex, b_pose,
                                      rotation_b_to_a, translation_b_in_a,
                                      request, contact_meta, result);
    const bool is_overlap_test9 = result.isCollision();

    EXPECT_EQ(is_overlap_test1, is_overlap);
    EXPECT_EQ(is_overlap_test2, is_overlap);
    EXPECT_EQ(is_overlap_test3, is_overlap);
    EXPECT_EQ(is_overlap_test4, is_overlap);
    EXPECT_EQ(is_overlap_test5, is_overlap);
    EXPECT_EQ(is_overlap_test6, is_overlap);
    EXPECT_EQ(is_overlap_test7, is_overlap);
    EXPECT_EQ(is_overlap_test8, is_overlap);
    EXPECT_EQ(is_overlap_test9, is_overlap);
  }
}

template <typename S>
void test_tetrahedron_box() {
  const int test_n = 10;
  for (int i = 0; i < test_n; i++) {
    const auto tetrahedron = test::generateRandomTetrahedron<S>();
    const auto box = test::generateRandomBox<S>();
    testTetrahedronWithBox<S>(tetrahedron, box);
  }
}

template <typename S>
void test_tetrahedron_triangle() {
  const int test_n = 10;
  for (int i = 0; i < test_n; i++) {
    const auto tetrahedron = test::generateRandomTetrahedron<S>();
    const auto triangle = test::generateRandomTriangleP<S>();
    testTetrahedronWithTriangleP<S>(tetrahedron, triangle);
  }
}

template <typename S>
void test_tetrahedron_tetrahedron() {
  const int test_n = 10;
  for (int i = 0; i < test_n; i++) {
    const auto a = test::generateRandomTetrahedron<S>();
    const auto b = test::generateRandomTetrahedron<S>();
    testTetrahedronWithTetrahedron<S>(a, b);
  }
}

GTEST_TEST(FCL_TETRAHEDRON, collision_tetrahedron_box) {
  test_tetrahedron_box<float>();
  test_tetrahedron_box<double>();
}

GTEST_TEST(FCL_TETRAHEDRON, collision_tetrahedron_tetrahedron) {
  test_tetrahedron_tetrahedron<float>();
  test_tetrahedron_tetrahedron<double>();
}

GTEST_TEST(FCL_TETRAHEDRON, collision_tetrahedron_triangle) {
  test_tetrahedron_triangle<float>();
  test_tetrahedron_triangle<double>();
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
