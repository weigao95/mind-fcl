#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace heightmap {

std::shared_ptr<PointCloud> generateRandomPointCloud(std::size_t numPoints) {
  auto cloud = std::make_shared<fcl::octomap::Pointcloud>();
  cloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    cloud->push_back(point[0], point[1], point[2] + 1);
  }
  return cloud;
}

std::shared_ptr<PointCloud> generateFullPointCloud() {
  auto cloud = std::make_shared<PointCloud>();

  // insert some points
  for (int x = -100; x < 100; x++) {
    for (int y = -100; y < 100; y++) {
      for (int z = 0; z < 20; z++) {
        cloud->push_back(octomap::point3d(x * 0.001, y * 0.001, z * 0.05));
      }
    }
  }

  return cloud;
}

template <typename S, typename Shape>
void heightmapCollisionRandomTestByBottomLayer(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  // Construct the solver
  const fcl::CollisionRequest<S> request{UINT_MAX};
  fcl::detail::GJKSolver<S> narrowphase_solver;
  narrowphase_solver.gjk_tolerance = request.binaryCollisionTolerance();
  narrowphase_solver.epa_tolerance = request.distanceTolerance();
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);

  // Check the collision between the shape and the small boxes that
  // make up the heightmap
  {
    fcl::CollisionResult<S> result;
    solver.HeightMapShapeIntersectByBottomLayer(&hm_geometry, shape, tf1, tf2,
                                                request, result);
    const auto& contacts = result.getContacts();
    std::set<int> contact_set;
    for (const auto& contact : contacts) {
      contact_set.insert(contact.b1);
    }

    {
      // Try with traverse
      fcl::CollisionResult<S> result_traverse;
      solver.HeightMapShapeIntersectByTraversal(&hm_geometry, shape, tf1, tf2,
                                                request, result_traverse);
      EXPECT_EQ(result_traverse.numContacts(), result.numContacts());
    }

    {
      // Try with fcl::collide
      fcl::CollisionResult<S> result_fcl;
      fcl::collide(&hm_geometry, tf1, &shape, tf2, request, result_fcl);
      EXPECT_EQ(result_fcl.numContacts(), result.numContacts());
    }

    // std::cerr << "The # of contacts in flat heightmap test "
    //           << contact_set.size() << std::endl;

    const FlatHeightMap<S>& bottom_map = hm_geometry.raw_heightmap()->bottom();
    for (uint16_t y = 0; y < bottom_map.full_shape_y(); y++) {
      for (uint16_t x = 0; x < bottom_map.full_shape_x(); x++) {
        Box<S> box;
        Transform3<S> box_tf;
        const heightmap::Pixel pixel = {x, y};
        bool ok = bottom_map.pixelToBox(pixel, box, box_tf);
        if (!ok) continue;
        box.computeLocalAABB();
        bool is_collision_with_bin = narrowphase_solver.shapeIntersect(
            box, tf1 * box_tf, shape, tf2, nullptr);

        const int code = heightmap::encodePixel(pixel);
        bool has_this_contact = (contact_set.find(code) != contact_set.end());
        EXPECT_EQ(has_this_contact, is_collision_with_bin);
      }
    }
  }
}

template <typename S, typename Shape>
void heightmapCollisionRandomTestByTraversal(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  // Construct the solver
  const fcl::detail::GJKSolver<S> narrowphase_solver;
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);

  // Check the collision between the shape and the small boxes that
  // make up the heightmap
  {
    const fcl::CollisionRequest<S> request{UINT_MAX};
    fcl::CollisionResult<S> result;
    solver.HeightMapShapeIntersectByTraversal(&hm_geometry, shape, tf1, tf2,
                                              request, result);
    const auto& contacts = result.getContacts();
    std::set<int> contact_set;
    for (const auto& contact : contacts) {
      contact_set.insert(contact.b1);
    }

    // std::cerr << "The # of contacts in flat heightmap test "
    //           << contact_set.size() << std::endl;

    const FlatHeightMap<S>& bottom_map = hm_geometry.raw_heightmap()->bottom();
    for (uint16_t y = 0; y < bottom_map.full_shape_y(); y++) {
      for (uint16_t x = 0; x < bottom_map.full_shape_x(); x++) {
        Box<S> box;
        Transform3<S> box_tf;
        const heightmap::Pixel pixel = {x, y};
        bool ok = bottom_map.pixelToBox(pixel, box, box_tf);
        if (!ok) continue;
        box.computeLocalAABB();
        bool is_collision_with_bin = narrowphase_solver.shapeIntersect(
            box, tf1 * box_tf, shape, tf2, nullptr);

        const int code = heightmap::encodePixel(pixel);
        bool has_this_contact = (contact_set.find(code) != contact_set.end());
        EXPECT_EQ(has_this_contact, is_collision_with_bin);
      }
    }
  }
}

template <typename S>
void heightmapCollisionRandomTests() {
  std::vector<fcl::HeightMapCollisionGeometry<S>> hm_geometries;
  {
    const auto points = generateRandomPointCloud(100000);
    auto heightMap = std::make_shared<LayeredHeightMap<S>>(0.002, 512);
    heightMap->updateHeightsByPointCloud3D(*points);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    const auto points = generateRandomPointCloud(1000000);
    auto heightMap =
        std::make_shared<LayeredHeightMap<S>>(0.002, 0.004, 256, 512);
    heightMap->updateHeightsByPointCloud3D(*points);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    const auto points = generateRandomPointCloud(500000);
    auto heightMap =
        std::make_shared<LayeredHeightMap<S>>(0.001, 0.002, 512, 512);
    heightMap->updateHeightsByPointCloud3D(*points);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  {
    const auto points = generateRandomPointCloud(300000);
    auto heightMap =
        std::make_shared<LayeredHeightMap<S>>(0.002, 0.001, 256, 512);
    heightMap->updateHeightsByPointCloud3D(*points);
    auto geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
    geometry.computeLocalAABB();
    hm_geometries.push_back(geometry);
  }

  // Test loop
  std::array<S, 6> extent{-1, -1, -1, 1, 1, 1};
  for (const auto& geometry : hm_geometries) {
    for (int i = 0; i < 10; i++) {
      Eigen::Transform<S, 3, Eigen::Isometry> tf1;
      test::generateRandomTransform(extent, tf1);
      Eigen::Transform<S, 3, Eigen::Isometry> tf2;
      test::generateRandomTransform(extent, tf2);

      fcl::Box<S> box(0.4, 0.4, 0.4);
      box.computeLocalAABB();
      heightmapCollisionRandomTestByBottomLayer<S>(geometry, box, tf1, tf2);
      heightmapCollisionRandomTestByTraversal<S>(geometry, box, tf1, tf2);

      fcl::Sphere<S> sphere(0.4);
      sphere.computeLocalAABB();
      heightmapCollisionRandomTestByBottomLayer<S>(geometry, sphere, tf1, tf2);
      heightmapCollisionRandomTestByTraversal<S>(geometry, sphere, tf1, tf2);

      fcl::Cylinder<S> cylinder(0.4, 0.5);
      cylinder.computeLocalAABB();
      heightmapCollisionRandomTestByBottomLayer<S>(geometry, cylinder, tf1,
                                                   tf2);
      heightmapCollisionRandomTestByTraversal<S>(geometry, cylinder, tf1, tf2);
    }
  }
}

template <typename S, typename Shape>
fcl::CollisionResult<S> runHeightMapCollision(
    const fcl::HeightMapCollisionGeometry<S>& map, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2, bool traverse = false) {
  /// Check the collision between the heightmap and the shape
  fcl::detail::GJKSolver<S> narrowphase_solver;
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);
  const fcl::CollisionRequest<S> request{UINT_MAX};
  fcl::CollisionResult<S> result;
  if (traverse) {
    solver.ShapeHeightMapIntersectByTraversal(shape, &map, tf2, tf1, request,
                                              result);
  } else {
    solver.HeightMapShapeIntersectByBottomLayer(&map, shape, tf1, tf2, request,
                                                result);
  }
  return result;
}

template <typename S>
void heightmapFullCloudCollisionTests() {
  const auto points = generateFullPointCloud();
  float resolution = 0.002;
  auto heightMap = std::make_shared<LayeredHeightMap<S>>(resolution, 256);
  heightMap->updateHeightsByPointCloud3D(*points);
  auto g = fcl::HeightMapCollisionGeometry<S>(heightMap);
  g.computeLocalAABB();

  fcl::Box<S> box(0.02, 0.02, 0.02);
  box.computeLocalAABB();
  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{1, 1, 0.5};
    auto result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, false);
    EXPECT_FALSE(result.isCollision());
    result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, true);
    EXPECT_FALSE(result.isCollision());
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.11 + 1e-6, 0.11 + 1e-6, 0.5};
    auto result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, false);
    EXPECT_FALSE(result.isCollision());
    result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, true);
    EXPECT_FALSE(result.isCollision());
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.1 + 1e-6, 0.1 + 1e-6, 0};
    auto result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 25u);
    result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 25u);
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.09 + 1e-6, 0.09 + 1e-6, 0};
    auto result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 100u);
    result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 100u);
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.0, 0.0, 0.0};
    auto result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 100u);
    result = runHeightMapCollision<S>(
        g, box, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf, true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 100u);
  }

  fcl::Sphere<S> sphere(0.01);
  sphere.computeLocalAABB();
  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{1, 1, 0.5};
    auto result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        false);
    EXPECT_FALSE(result.isCollision());
    result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        true);
    EXPECT_FALSE(result.isCollision());
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.11, 0.11, 0.5};
    auto result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        false);
    EXPECT_FALSE(result.isCollision());
    result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        true);
    EXPECT_FALSE(result.isCollision());
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.1 + 1e-6, 0.1 + 1e-6, 0};
    auto result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 22u);
    result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 22u);
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.09 + 1e-6, 0.09 + 1e-6, 0};
    auto result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 92u);
    result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 92u);
  }

  {
    Eigen::Transform<S, 3, Eigen::Isometry> tf;
    tf.linear().setIdentity();
    tf.translation() = Vector3<S>{0.0, 0.0, 0.0};
    auto result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        false);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 88u);
    result = runHeightMapCollision<S>(
        g, sphere, Eigen::Transform<S, 3, Eigen::Isometry>::Identity(), tf,
        true);
    EXPECT_TRUE(result.isCollision());
    EXPECT_EQ(result.numContacts(), 88u);
  }
}

}  // namespace heightmap
}  // namespace fcl

GTEST_TEST(HeightMapShapeCollision, RandomTest) {
  fcl::heightmap::heightmapCollisionRandomTests<float>();
  fcl::heightmap::heightmapCollisionRandomTests<double>();
}

GTEST_TEST(HeightMapShapeCollision, FullCloudTest) {
  fcl::heightmap::heightmapFullCloudCollisionTests<float>();
  fcl::heightmap::heightmapFullCloudCollisionTests<double>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
