//
// Created by Wei Gao on 2023/3/18.
//
#include <gtest/gtest.h>

#include "fcl/geometry/shape/shape_swept_volume.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/narrowphase/collision.h"

namespace fcl {

template <typename S>
void testSphereSweptVolume() {
  const S radius = 0.1;
  const S displacement_z = 0.6;
  auto sphere = std::make_shared<Sphere<S>>(radius);
  sphere->computeLocalAABB();
  Vector3<S> displacement;
  displacement.setZero();
  displacement.z() = displacement_z;
  ShapeTranslationSweptVolume<S> swept_volume =
      ShapeTranslationSweptVolume<S>::template Make<fcl::Sphere<S>>(
          sphere, displacement);
  EXPECT_NEAR(swept_volume.movement().z(), displacement.z(), 1e-6);

  std::vector<Vector3<S>> bounded_vertices =
      swept_volume.getBoundVertices(Transform3<S>::Identity());
  EXPECT_EQ(bounded_vertices.size(), 24u);

  Transform3<S> tf1;
  tf1.setIdentity();
  Transform3<S> tf2;
  tf2.setIdentity();
  CollisionRequest<S> request;
  CollisionResult<S> result;
  fcl::collide(&swept_volume, tf1, sphere.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());

  tf2.translation().x() = 2 * radius + 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, sphere.get(), tf2, request, result);
  EXPECT_FALSE(result.isCollision());

  tf2.translation().x() = 2 * radius - 1e-3;
  tf2.translation().z() = displacement_z * 0.5;
  request.useDefaultPenetration();
  result.clear();
  fcl::collide(&swept_volume, tf1, sphere.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());
  EXPECT_NEAR(1e-3, result.getContact(0).penetration_depth, 1e-5);

  tf2.setIdentity();
  tf2.translation().z() = displacement_z + 2 * radius + 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, sphere.get(), tf2, request, result);
  EXPECT_FALSE(result.isCollision());

  tf2.setIdentity();
  tf2.translation().z() = displacement_z + 2 * radius - 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, sphere.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());
  EXPECT_NEAR(1e-3, result.getContact(0).penetration_depth, 1e-5);
}

template <typename S>
void testBoxSweptVolume() {
  const S box_halfsize = 0.1;
  const S box_size = 2 * box_halfsize;
  const S displacement_z = 0.6;
  auto box = std::make_shared<Box<S>>(box_size, box_size, box_size);
  box->computeLocalAABB();

  // Make the swept volume
  Vector3<S> displacement;
  displacement.setZero();
  displacement.z() = displacement_z;
  ShapeTranslationSweptVolume<S> swept_volume =
      ShapeTranslationSweptVolume<S>::template Make<fcl::Box<S>>(box,
                                                                 displacement);
  EXPECT_NEAR(swept_volume.movement().z(), displacement.z(), 1e-6);

  Transform3<S> tf1;
  tf1.setIdentity();
  Transform3<S> tf2;
  tf2.setIdentity();
  CollisionRequest<S> request;
  CollisionResult<S> result;
  fcl::collide(&swept_volume, tf1, box.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());

  tf2.translation().x() = 2 * box_halfsize + 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, box.get(), tf2, request, result);
  EXPECT_FALSE(result.isCollision());

  tf2.translation().x() = 2 * box_halfsize - 1e-3;
  tf2.translation().z() = displacement_z * 0.5;
  request.useDefaultPenetration();
  result.clear();
  fcl::collide(&swept_volume, tf1, box.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());
  EXPECT_NEAR(1e-3, result.getContact(0).penetration_depth, 1e-5);

  tf2.setIdentity();
  tf2.translation().z() = displacement_z + 2 * box_halfsize + 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, box.get(), tf2, request, result);
  EXPECT_FALSE(result.isCollision());

  tf2.setIdentity();
  tf2.translation().z() = displacement_z + 2 * box_halfsize - 1e-3;
  result.clear();
  fcl::collide(&swept_volume, tf1, box.get(), tf2, request, result);
  EXPECT_TRUE(result.isCollision());
  EXPECT_NEAR(1e-3, result.getContact(0).penetration_depth, 1e-5);
}

}  // namespace fcl

GTEST_TEST(SweptVolumeGeometryTest, SphereTest) {
  fcl::testSphereSweptVolume<float>();
  fcl::testSphereSweptVolume<double>();
}

GTEST_TEST(SweptVolumeGeometryTest, BoxTest) {
  fcl::testBoxSweptVolume<float>();
  fcl::testBoxSweptVolume<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
