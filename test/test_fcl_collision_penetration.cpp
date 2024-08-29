//
// Created by mech-mind_gw on 3/21/2023.
//
#include <gtest/gtest.h>

#include <iostream>

#include "create_primitive_mesh.h"
#include "fcl/fcl.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename Shape1, typename Shape2>
void testCollisionShapePair(const Shape1& shape_1, const Shape2& shape_2,
                            const Transform3<typename Shape1::S>& tf_1,
                            const Transform3<typename Shape1::S>& tf_2) {
  static_assert(std::is_same<typename Shape1::S, typename Shape2::S>::value,
                "Shape must be declared with same scalar type");
  using S = typename Shape1::S;
  CollisionRequest<S> no_pd_request;
  no_pd_request.disablePenetration();
  CollisionResult<S> no_pd_result;
  fcl::collide(&shape_1, tf_1, &shape_2, tf_2, no_pd_request, no_pd_result);
  if (!no_pd_result.isCollision()) return;

  // Test conditions
  CollisionRequest<S> request_default_pd;
  request_default_pd.useDefaultPenetration();
  CollisionRequest<S> request_unit_z_pd;
  request_unit_z_pd.useDirectedPenetration(Vector3<S>::UnitZ());
  CollisionRequest<S> request_incremental_unit_z_pd;
  request_incremental_unit_z_pd.useIncrementalMinimumDistancePenetration(
      Vector3<S>::UnitZ());

  // Start test loop
  for (const auto& request_pd :
       {request_default_pd, request_unit_z_pd, request_incremental_unit_z_pd}) {
    CollisionResult<S> result_pd;
    fcl::collide(&shape_1, tf_1, &shape_2, tf_2, request_pd, result_pd);
    EXPECT_TRUE(result_pd.isCollision());

    // Now there is a collision
    for (std::size_t i = 0; i < result_pd.numContacts(); i++) {
      const Contact<S>& contact_i = result_pd.getContact(i);
      if (!(contact_i.penetration_depth >= S(0.0))) {
        std::cout << tf_1.matrix() << std::endl;
        std::cout << "Now is tf2" << std::endl;
        std::cout << tf_2.matrix() << std::endl;
        CollisionResult<S> result_debug;
        fcl::collide(&shape_1, tf_1, &shape_2, tf_2, request_pd, result_debug);
      }
      EXPECT_TRUE(contact_i.penetration_depth >= S(0.0));
      const Vector3<S> shape_2_escape =
          (contact_i.penetration_depth * S(1.05)) * contact_i.normal;
      EXPECT_NEAR(
          (shape_2_escape - contact_i.shape2_escape_movement(1.05)).norm(), 0,
          S(1e-6));

      // Move shape2
      {
        Transform3<S> tf_2_moved = tf_2;
        tf_2_moved.translation() += shape_2_escape;
        fcl::CollisionResult<S> moved_result;
        fcl::collide(&shape_1, tf_1, &shape_2, tf_2_moved, no_pd_request,
                     moved_result);
        EXPECT_TRUE(!moved_result.isCollision());
      }

      // Move shape1
      {
        Transform3<S> tf_1_moved = tf_1;
        tf_1_moved.translation() += contact_i.shape1_escape_movement(1.05);
        fcl::CollisionResult<S> moved_result;
        fcl::collide(&shape_1, tf_1_moved, &shape_2, tf_2, no_pd_request,
                     moved_result);
        EXPECT_TRUE(!moved_result.isCollision());
      }
    }
    // Done with all contact
  }
}

template <typename Shape>
void sameTypeShapePairTest(const Shape& shape, std::size_t test_n) {
  using S = typename Shape::S;
  std::vector<Transform3<S>> tf1_vec, tf2_vec;
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  test::generateRandomTransformVector(extent, test_n, tf1_vec);
  test::generateRandomTransformVector(extent, test_n, tf2_vec);
  for (std::size_t i = 0; i < test_n; i++) {
    testCollisionShapePair(shape, shape, tf1_vec[i], tf2_vec[i]);
  }
}

template <typename Shape1, typename Shape2>
void differentTypeShapePairTest(const Shape1& shape1, const Shape2& shape2,
                                std::size_t test_n) {
  using S = typename Shape1::S;
  std::vector<Transform3<S>> tf1_vec, tf2_vec;
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  test::generateRandomTransformVector(extent, test_n, tf1_vec);
  test::generateRandomTransformVector(extent, test_n, tf2_vec);
  for (std::size_t i = 0; i < test_n; i++) {
    testCollisionShapePair<Shape1, Shape2>(shape1, shape2, tf1_vec[i],
                                           tf2_vec[i]);
    testCollisionShapePair<Shape2, Shape1>(shape2, shape1, tf1_vec[i],
                                           tf2_vec[i]);
  }
}

template <typename Shape1, typename Shape2>
void differentTypeShapePairTestInstance(const Shape1& shape1,
                                        const Shape2& shape2) {
  using S = typename Shape1::S;
  Transform3<S> shape1_pose, shape2_pose;
  {
    shape1_pose.setIdentity();
    shape1_pose.translation().x() = 0.1957880244590342045;
    shape1_pose.translation().y() = 0.2207083962857722681;
    shape1_pose.translation().z() = 0.1890375458635389916;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        -0.05616368847400071002, 0.02208880723252211811, 0.9981772010480102209);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        0.8777983833572846617, 0.477454328244468984, 0.03882476807349017484);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        -0.4757264321774992144, 0.8783788755631976031, -0.04620511543108424268);
    shape1_pose.linear().matrix() = rotation_matrix;
  }

  {
    shape2_pose.setIdentity();
    shape2_pose.translation().x() = 0.2755832296796142944;
    shape2_pose.translation().y() = 0.2450562855228781145;
    shape2_pose.translation().z() = 0.1936604700051248185;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.2156778687337976008, 0.3274077108955088522, 0.9199387195811521423);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        -0.7799624058762353318, -0.5090555864333072833, 0.3640344150500010034);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        0.5874876188977956604, -0.7960317837669577123, 0.1455736819435227425);
    shape2_pose.linear().matrix() = rotation_matrix;
  }
  testCollisionShapePair<Shape1, Shape2>(shape1, shape2, shape1_pose,
                                         shape2_pose);
}

// Pairs of same shape
template <typename S>
void spherePairTest(std::size_t test_n) {
  const S radius = 0.4;
  Sphere<S> shape(radius);
  fcl::sameTypeShapePairTest<Sphere<S>>(shape, test_n);
}

template <typename S>
void boxPairTest(std::size_t test_n) {
  const S half_size = 0.4;
  Box<S> shape(2 * half_size, 2 * half_size, 2 * half_size);
  fcl::sameTypeShapePairTest<Box<S>>(shape, test_n);
}

template <typename S>
void capsulePairTest(std::size_t test_n) {
  const S radius = 0.4;
  const S lz = 0.8;
  Capsule<S> shape(radius, lz);
  fcl::sameTypeShapePairTest<Capsule<S>>(shape, test_n);
}

template <typename S>
void cylinderPairTest(std::size_t test_n) {
  const S radius = 0.4;
  const S lz = 0.8;
  Cylinder<S> shape(radius, lz);
  fcl::sameTypeShapePairTest<Cylinder<S>>(shape, test_n);
}

template <typename S>
void convexPairTest(std::size_t test_n) {
  const S radius = 0.4;
  auto convex = std::make_shared<Convex<S>>(test::makeSphereConvex<S>(radius));
  convex->computeLocalAABB();
  fcl::sameTypeShapePairTest<Convex<S>>(*convex, test_n);
}

// Pairs of different shape
template <typename S>
void sphereBoxTest(std::size_t test_n) {
  const S radius = 0.4;
  Sphere<S> shape1(radius);
  const S half_size = 0.4;
  Box<S> shape2(2 * half_size, 2 * half_size, 2 * half_size);
  fcl::differentTypeShapePairTest(shape1, shape2, test_n);
}

template <typename S>
void sphereCapsuleTest(std::size_t test_n) {
  const S radius = 0.4;
  Sphere<S> shape1(radius);

  const S lz = 0.8;
  Capsule<S> shape2(radius, lz);
  // fcl::differentTypeShapePairTestInstance(shape1, shape2);
  fcl::differentTypeShapePairTest(shape1, shape2, test_n);
}

template <typename S>
void sphereCylinderTest(std::size_t test_n) {
  const S radius = 0.4;
  Sphere<S> shape1(radius);

  const S lz = 0.8;
  Cylinder<S> shape2(radius, lz);
  fcl::differentTypeShapePairTest(shape1, shape2, test_n);
}

}  // namespace fcl

// Same shape pairs
GTEST_TEST(CollisionPenetrationTest, SpherePairTest) {
  std::size_t test_n = 1000u;
  fcl::spherePairTest<float>(test_n);
  fcl::spherePairTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, BoxPairTest) {
  std::size_t test_n = 1000u;
  fcl::boxPairTest<float>(test_n);
  fcl::boxPairTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, CylinderPairTest) {
  std::size_t test_n = 1000u;
  fcl::cylinderPairTest<float>(test_n);
  fcl::cylinderPairTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, CapsulePairTest) {
  std::size_t test_n = 1000u;
  fcl::capsulePairTest<float>(test_n);
  fcl::capsulePairTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, ConvexPairTest) {
  std::size_t test_n = 1000u;
  fcl::convexPairTest<float>(test_n);
  fcl::convexPairTest<double>(test_n);
}

// Different shape pairs
GTEST_TEST(CollisionPenetrationTest, SphereBoxTest) {
  std::size_t test_n = 1000u;
  fcl::sphereBoxTest<float>(test_n);
  fcl::sphereBoxTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, SphereCapsuleTest) {
  std::size_t test_n = 1000u;
  fcl::sphereCapsuleTest<float>(test_n);
  fcl::sphereCapsuleTest<double>(test_n);
}

GTEST_TEST(CollisionPenetrationTest, SphereCylinderTest) {
  std::size_t test_n = 1000u;
  fcl::sphereCylinderTest<float>(test_n);
  fcl::sphereCylinderTest<double>(test_n);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}