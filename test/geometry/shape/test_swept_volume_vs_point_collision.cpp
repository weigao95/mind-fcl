//
// Created by mech-mind_gw on 3/20/2023.
//
#include <gtest/gtest.h>

#include "create_primitive_mesh.h"
#include "fcl/geometry/shape/shape_swept_volume.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/translational_ccd.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename Shape1, typename Shape2>
bool shapeSweptTest(const std::shared_ptr<const Shape1>& shape1,
                    const std::shared_ptr<const Shape2>& shape2,
                    const Transform3<typename Shape1::S>& tf1,
                    const Transform3<typename Shape1::S>& tf2,
                    const Vector3<typename Shape1::S>& shape1_sweep_in_world,
                    std::size_t discrete_points_size) {
  using S = typename Shape1::S;
  Transform3<S> world2shape1 = tf1.inverse(Eigen::Isometry);
  Vector3<S> shape1_sweep_in_shape1 =
      world2shape1.linear() * shape1_sweep_in_world;
  auto swept_volume = std::make_shared<ShapeTranslationSweptVolume<S>>(
      ShapeTranslationSweptVolume<S>::template Make<Shape1>(
          shape1, shape1_sweep_in_shape1));
  swept_volume->computeLocalAABB();

  // Compute the swept volume collision
  CollisionRequest<S> request;
  CollisionResult<S> result;
  fcl::collide(swept_volume.get(), tf1, shape2.get(), tf2, request, result);
  const bool swept_collide = result.isCollision();

  bool collide_by_point = false;
  for (std::size_t i = 0; i < discrete_points_size; i++) {
    const S t = S(i) / S(discrete_points_size - 1);
    Vector3<S> displacement_t = t * shape1_sweep_in_world;
    Transform3<S> tf1_t = tf1;
    tf1_t.translation() += displacement_t;

    result.clear();
    fcl::collide(shape1.get(), tf1_t, shape2.get(), tf2, request, result);
    if (result.isCollision()) {
      collide_by_point = true;
      break;
    }
  }

  // Test by fcl
  CollisionObject<S> o1(
      shape1, tf1, typename CollisionObject<S>::GeometryLocalAABBComputed());
  CollisionObject<S> o2(
      shape2, tf2, typename CollisionObject<S>::GeometryLocalAABBComputed());
  result.clear();
  translationalContinuousCollision(&o1, &o2, shape1_sweep_in_world,
                                   request, result);
  const bool collide_by_fcl_interface = result.isCollision();
  EXPECT_EQ(collide_by_fcl_interface, swept_collide);

  // If colliding, then their AABB must intersect
  if (swept_collide) {
    const AABB<S>& swept_aabb = swept_volume->aabb_local;
    const AABB<S>& shape_aabb = shape2->aabb_local;
    OBB<S> swept_obb, shape_obb;
    convertBV(swept_aabb, tf1, swept_obb);
    convertBV(shape_aabb, tf2, shape_obb);
    EXPECT_TRUE(swept_obb.overlap(shape_obb));
  }

  // Done
  return (collide_by_point == swept_collide);
}

template <typename Shape1, typename Shape2>
std::size_t compareCollisionResultWithPointCollide(
    const std::shared_ptr<const Shape1>& shape1,
    const std::shared_ptr<const Shape2>& shape2,
    const std::vector<Transform3<typename Shape1::S>>& tf1_vec,
    const std::vector<Transform3<typename Shape1::S>>& tf2_vec,
    const std::vector<Vector3<typename Shape1::S>>& shape1_sweep_vec) {
  EXPECT_EQ(tf1_vec.size(), tf2_vec.size());
  std::size_t n_mismatch = 0;
  for (std::size_t i = 0; i < tf1_vec.size(); i++) {
    bool result_match = shapeSweptTest(shape1, shape2, tf1_vec[i], tf2_vec[i],
                                       shape1_sweep_vec[i], 20);
    if (!result_match) {
      n_mismatch += 1;
    }
  }
  return n_mismatch;
}

template <typename S>
void testSweptVolumeWithPointCollide() {
  // Generate the collision of random pose
  const std::size_t test_n = 100;
  std::vector<Transform3<S>> tf1_vec, tf2_vec;
  std::vector<Vector3<S>> shape1_displacement;
  std::array<S, 6> extent{-1, -1, -1, 1, 1, 1};
  for (std::size_t test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> shape1_pose, shape2_pose, displacement;
    test::generateRandomTransform(extent, shape1_pose);
    test::generateRandomTransform(extent, shape2_pose);
    test::generateRandomTransform(extent, displacement);
    tf1_vec.emplace_back(shape1_pose);
    tf2_vec.emplace_back(shape2_pose);
    shape1_displacement.emplace_back(S(0.3) * displacement.translation());
  }

  const S radius = 0.4;
  const S lz = 0.8;
  auto capsule = std::make_shared<Capsule<S>>(radius, lz);
  capsule->computeLocalAABB();

  auto sphere = std::make_shared<Sphere<S>>(radius);
  sphere->computeLocalAABB();

  auto cylinder = std::make_shared<Cylinder<S>>(radius, lz);
  cylinder->computeLocalAABB();

  auto box = std::make_shared<Box<S>>(0.4, 0.6, 0.8);
  box->computeLocalAABB();

  auto convex = std::make_shared<Convex<S>>(test::makeSphereConvex<S>(radius));
  convex->computeLocalAABB();

  // Sphere vs sphere
  {
    auto n_mismatch =
        compareCollisionResultWithPointCollide<Sphere<S>, Sphere<S>>(
            sphere, sphere, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Box vs Box
  {
    auto n_mismatch = compareCollisionResultWithPointCollide<Box<S>, Box<S>>(
        box, box, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Cylinder vs capsule
  {
    auto n_mismatch =
        compareCollisionResultWithPointCollide<Cylinder<S>, Capsule<S>>(
            cylinder, capsule, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Convex vs capsule
  {
    auto n_mismatch =
        compareCollisionResultWithPointCollide<Convex<S>, Capsule<S>>(
            convex, capsule, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Convex vs box
  {
    auto n_mismatch = compareCollisionResultWithPointCollide<Convex<S>, Box<S>>(
        convex, box, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Convex vs convex
  {
    auto n_mismatch =
        compareCollisionResultWithPointCollide<Convex<S>, Convex<S>>(
            convex, convex, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }

  // Convex vs octree
  {
    const S resolution = 0.001;  // 1 mm
    auto octree = std::shared_ptr<const octomap::OcTree>(
        test::generateOcTree(resolution));
    auto fcl_octree = std::make_shared<OcTree<S>>(octree);
    fcl_octree->computeLocalAABB();

    auto n_mismatch =
        compareCollisionResultWithPointCollide<Convex<S>, OcTree<S>>(
            convex, fcl_octree, tf1_vec, tf2_vec, shape1_displacement);
    EXPECT_TRUE(n_mismatch == 0);
  }
}

}  // namespace fcl

GTEST_TEST(SweptVolumePointCollisionTest, PointCollisionTest) {
  fcl::testSweptVolumeWithPointCollide<float>();
  fcl::testSweptVolumeWithPointCollide<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}