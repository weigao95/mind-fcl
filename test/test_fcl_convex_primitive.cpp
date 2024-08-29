//
// Created by wei on 2/1/22.
//
#include "create_primitive_mesh.h"
#include "test_fcl_utility.h"
#include <gtest/gtest.h>

template <typename T>
void testConvexSupport(const fcl::Convex<T>& convex) {
  const int test_n = 1000;
  const auto& convex_vertices = convex.getVertices();
  for(auto i = 0; i < test_n; i++) {
    fcl::Vector3<T> direction;
    direction.setRandom();
    if (direction.norm() < 1e-6)
      continue;

    // Unit direction
    direction.normalize();
    const auto& vertex_neighbor = convex.findExtremeVertex(direction);

    // The naive vertex
    fcl::Vector3<T> vertex_naive = convex_vertices[0];
    T max_value = direction.dot(vertex_naive);
    for(std::size_t j = 0; j < convex_vertices.size(); j++) {
      auto dot_j = convex_vertices[j].dot(direction);
      if (dot_j > max_value) {
        max_value = dot_j;
        vertex_naive = convex_vertices[j];
      }
    }

    // The result should match
    EXPECT_TRUE((vertex_naive - vertex_neighbor).norm() < 1e-7);
  }
}

template <typename T>
void sphereConvexTestInstance() {
  T radius = 0.1;
  auto sphere_1 = std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));
  auto sphere_2 = std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));

  // Test on given input
  fcl::Matrix3<T> rotation_1, rotation_2;
  fcl::Vector3<T> translation_1, translation_2;
  rotation_1.row(0) = fcl::Vector3<T>(-0.649263, 0.203951, -0.732708);
  rotation_1.row(1) = fcl::Vector3<T>(0.678622, 0.590328, -0.437018);
  rotation_1.row(2) = fcl::Vector3<T>(0.343408, -0.780972, -0.521684);
  translation_1 = fcl::Vector3<T>(0.494675, -0.00619526, -0.239784);

  rotation_2.row(0) = fcl::Vector3<T>(0.413224, 0.822318, -0.391202);
  rotation_2.row(1) = fcl::Vector3<T>(-0.783547, 0.102179, -0.612873);
  rotation_2.row(2) = fcl::Vector3<T>(-0.464004, 0.559778, 0.686548);
  translation_2 = fcl::Vector3<T>(0.570727, 0.105625, -0.285806);

  fcl::Transform3<T> body_1_pose, body_2_pose;
  body_1_pose.linear() = rotation_1;
  body_1_pose.translation() = translation_1;
  body_2_pose.linear() = rotation_2;
  body_2_pose.translation() = translation_2;

  // Make the object
  fcl::CollisionObject<T> object_1(sphere_1, body_1_pose);
  fcl::CollisionObject<T> object_2(sphere_2, body_2_pose);
  auto distance =
      (body_2_pose.translation() - body_1_pose.translation()).norm();
  bool should_collide = distance < 2 * radius;

  // Run detection
  fcl::CollisionRequest<T> request;
  fcl::CollisionResult<T> result;
  fcl::collide(&object_1, &object_2, request, result);
  EXPECT_TRUE(result.isCollision() == should_collide);
}

template <typename T>
void sphereTestRandomPose() {
  T radius = 0.1;
  auto sphere_1 = std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));
  auto sphere_2 = std::make_shared<fcl::Convex<T>>(fcl::test::makeSphereConvex<T>(radius));

  // The extent used to generate the pose
  T extent = radius * 2;
  std::array<T, 6> xyz_extent{
      -extent, -extent, -extent, extent, extent, extent};

  int test_n = 1e3;
  for(auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<T> body_1_pose, body_2_pose;
    fcl::test::generateRandomTransform(xyz_extent, body_1_pose);
    fcl::test::generateRandomTransform(xyz_extent, body_2_pose);
    auto distance =
        (body_2_pose.translation() - body_1_pose.translation()).norm();
    bool should_collide = distance <= 2 * radius;
    double tolerance = std::abs(distance - 2 * radius);

    // Make the object
    fcl::CollisionObject<T> object_1(sphere_1, body_1_pose);
    fcl::CollisionObject<T> object_2(sphere_2, body_2_pose);
    fcl::CollisionRequest<T> request;
    fcl::CollisionResult<T> result;
    fcl::collide(&object_1, &object_2, request, result);
    if(!result.isCollision() == should_collide && tolerance > 0.1 * radius) {
      std::cout << "Mismatch at tolerance " << tolerance << std::endl;
      std::cout << "Body 1 pose" << std::endl;
      std::cout << body_1_pose.matrix() << std::endl;
      std::cout << "Body 2 pose" << std::endl;
      std::cout << body_2_pose.matrix() << std::endl;
      std::cout << "******************************" << std::endl;
    }
  }
}

GTEST_TEST(FCL_Convex_Primitive, SupportTest) {
  double radius = 0.1;
  auto sphere_double = fcl::test::makeSphereConvex<double>(radius);
  auto sphere_float = fcl::test::makeSphereConvex<float>(radius);
  testConvexSupport(sphere_double);
  testConvexSupport(sphere_float);
}

GTEST_TEST(FCL_Convex_Primitive, SphereTest) {
  sphereConvexTestInstance<float>();
  sphereConvexTestInstance<double>();
}

GTEST_TEST(FCL_Convex_Primitive, SphereRandomTest) {
  sphereTestRandomPose<float>();
  sphereTestRandomPose<double>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
