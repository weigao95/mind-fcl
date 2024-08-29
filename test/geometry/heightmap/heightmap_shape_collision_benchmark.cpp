#include <chrono>
#include <vector>

#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace heightmap {

std::shared_ptr<const PointCloud> generateRandomPointCloud(
    std::size_t numPoints) {
  auto pointCloud = std::make_shared<PointCloud>();
  pointCloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    pointCloud->push_back(point[0], point[1], point[2] + 1);
  }
  return pointCloud;
}

template <typename S, typename Shape>
void heightmapCollisionRandomTest(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  const fcl::detail::GJKSolver<S> narrowphase_solver;
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);
  const fcl::CollisionRequest<S> request{UINT_MAX};
  fcl::CollisionResult<S> result;
  solver.HeightMapShapeIntersectByBottomLayer(&hm_geometry, shape, tf1, tf2,
                                              request, result);
}

template <typename S, typename Shape>
void heightmapCollisionRandomTestByTraversal(
    const fcl::HeightMapCollisionGeometry<S>& hm_geometry, const Shape& shape,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf1,
    const Eigen::Transform<S, 3, Eigen::Isometry>& tf2) {
  const fcl::detail::GJKSolver<S> narrowphase_solver;
  fcl::detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);
  const fcl::CollisionRequest<S> request{UINT_MAX};
  fcl::CollisionResult<S> result;
  solver.HeightMapShapeIntersectByTraversal(&hm_geometry, shape, tf1, tf2,
                                            request, result);
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

  std::array<S, 6> extent{-1, -1, -1, 1, 1, 1};
  std::vector<Eigen::Transform<S, 3, Eigen::Isometry>> v_tf1, v_tf2;
  for (int i = 0; i < 10; i++) {
    Eigen::Transform<S, 3, Eigen::Isometry> tf1;
    test::generateRandomTransform(extent, tf1);
    v_tf1.push_back(tf1);
    Eigen::Transform<S, 3, Eigen::Isometry> tf2;
    test::generateRandomTransform(extent, tf2);
    v_tf2.push_back(tf2);
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (const auto& geometry : hm_geometries) {
    for (const auto& tf1 : v_tf1) {
      for (const auto& tf2 : v_tf2) {
        fcl::Box<S> box(0.4, 0.4, 0.4);
        box.computeLocalAABB();
        heightmapCollisionRandomTest<S>(geometry, box, tf1, tf2);

        fcl::Sphere<S> sphere(0.4);
        sphere.computeLocalAABB();
        heightmapCollisionRandomTest<S>(geometry, sphere, tf1, tf2);

        fcl::Cylinder<S> cylinder(0.4, 0.5);
        cylinder.computeLocalAABB();
        heightmapCollisionRandomTest<S>(geometry, cylinder, tf1, tf2);
      }
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "HeightMapShapeIntersectByBottomLayer Time in ms: " << ms_time
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (const auto& geometry : hm_geometries) {
    for (const auto& tf1 : v_tf1) {
      for (const auto& tf2 : v_tf2) {
        fcl::Box<S> box(0.4, 0.4, 0.4);
        box.computeLocalAABB();
        heightmapCollisionRandomTestByTraversal<S>(geometry, box, tf1, tf2);

        fcl::Sphere<S> sphere(0.4);
        sphere.computeLocalAABB();
        heightmapCollisionRandomTestByTraversal<S>(geometry, sphere, tf1, tf2);

        fcl::Cylinder<S> cylinder(0.4, 0.5);
        cylinder.computeLocalAABB();
        heightmapCollisionRandomTestByTraversal<S>(geometry, cylinder, tf1,
                                                   tf2);
      }
    }
  }

  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "HeightMapShapeIntersectByTraversal Time in ms: " << ms_time
            << std::endl;
}

}  // namespace heightmap
}  // namespace fcl

//==============================================================================
int main() {
  std::cout << "Benchmark with float" << std::endl;
  fcl::heightmap::heightmapCollisionRandomTests<float>();
  std::cout << "Benchmark with double" << std::endl;
  fcl::heightmap::heightmapCollisionRandomTests<double>();
}
