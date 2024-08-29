//
// Created by mech-mind_gw on 3/20/2023.
//
#include <gtest/gtest.h>
#include "create_primitive_mesh.h"
#include "fcl/geometry/shape/shape_swept_volume.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename Shape>
bool compareCollisionResultOnPose(
    const Shape& shape2test,
    const ShapeTranslationSweptVolume<typename Shape::S>& sphere_z_swept,
    const Capsule<typename Shape::S>& capsule,
    const Transform3<typename Shape::S>& tf_shape,
    const Transform3<typename Shape::S>& tf_swept) {
  using S = typename Shape::S;
  Transform3<S> capsule2swept;
  capsule2swept.setIdentity();
  capsule2swept.translation().z() = S(0.5) * capsule.lz;
  Transform3<S> tf_capsule = tf_swept * capsule2swept;

  CollisionRequest<S> request;
  CollisionResult<S> result;
  fcl::collide(&shape2test, tf_shape, &capsule, tf_capsule, request, result);
  const bool capsule_collide = result.isCollision();

  result.clear();
  fcl::collide(&shape2test, tf_shape, &sphere_z_swept, tf_swept, request,
               result);
  const bool swept_collide = result.isCollision();
  EXPECT_EQ(swept_collide, capsule_collide);

  // If colliding, then their AABB must intersect
  if (capsule_collide) {
    const AABB<S>& swept_aabb = sphere_z_swept.aabb_local;
    const AABB<S>& shape_aabb = shape2test.aabb_local;
    OBB<S> swept_obb, shape_obb;
    convertBV(swept_aabb, tf_swept, swept_obb);
    convertBV(shape_aabb, tf_shape, shape_obb);
    EXPECT_TRUE(swept_obb.overlap(shape_obb));
  }

  return swept_collide == capsule_collide;
}

template <typename Shape>
void compareCollisionResult(
    const Shape& shape2test,
    const ShapeTranslationSweptVolume<typename Shape::S>& sphere_z_swept,
    const Capsule<typename Shape::S>& capsule,
    const std::vector<Transform3<typename Shape::S>>& tf_shape_vec,
    const std::vector<Transform3<typename Shape::S>>& tf_swept_vec) {
  EXPECT_EQ(tf_shape_vec.size(), tf_swept_vec.size());
  for (std::size_t i = 0; i < tf_shape_vec.size(); i++) {
    compareCollisionResultOnPose<Shape>(shape2test, sphere_z_swept, capsule,
                                        tf_shape_vec[i], tf_swept_vec[i]);
  }
}

template <typename S>
void testSphereSweptAsCapsule() {
  const S radius = 0.4;
  const S lz = 0.8;
  const std::size_t test_n = 1000;
  auto capsule = std::make_shared<Capsule<S>>(radius, lz);
  capsule->computeLocalAABB();

  auto sphere = std::make_shared<Sphere<S>>(radius);
  sphere->computeLocalAABB();

  Vector3<S> displacement;
  displacement.setZero();
  displacement.z() = lz;
  auto swept_volume = std::make_shared<ShapeTranslationSweptVolume<S>>(
      ShapeTranslationSweptVolume<S>::template Make<fcl::Sphere<S>>(
          sphere, displacement));
  swept_volume->computeLocalAABB();

  // Generate the collision of random pose
  std::vector<Transform3<S>> tf_shape, tf_swept;
  std::array<S, 6> extent{-1, -1, -1, 1, 1, 1};
  for (std::size_t test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> shape1_pose, shape2_pose;
    test::generateRandomTransform(extent, shape1_pose);
    test::generateRandomTransform(extent, shape2_pose);
    tf_shape.emplace_back(shape1_pose);
    tf_swept.emplace_back(shape2_pose);
  }

  // Test sphere
  compareCollisionResult<Sphere<S>>(*sphere, *swept_volume, *capsule, tf_shape,
                                    tf_swept);

  // Test capsule
  compareCollisionResult<Capsule<S>>(*capsule, *swept_volume, *capsule,
                                     tf_shape, tf_swept);

  // Test box
  Box<S> box2test(0.4, 0.6, 0.8);
  box2test.computeLocalAABB();
  compareCollisionResult<Box<S>>(box2test, *swept_volume, *capsule, tf_shape,
                                 tf_swept);

  // Test cylinder
  Cylinder<S> cylinder2test(radius, lz);
  cylinder2test.computeLocalAABB();
  compareCollisionResult<Cylinder<S>>(cylinder2test, *swept_volume, *capsule,
                                      tf_shape, tf_swept);

  // Test convex
  Convex<S> convex = fcl::test::makeSphereConvex<S>(radius);
  convex.computeLocalAABB();
  compareCollisionResult<Convex<S>>(convex, *swept_volume, *capsule, tf_shape,
                                    tf_swept);
}

}  // namespace fcl

GTEST_TEST(SweptVolumeSphereCapsuleTest, SphereSweptAsCapsuleTest) {
  fcl::testSphereSweptAsCapsule<float>();
  fcl::testSphereSweptAsCapsule<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
