#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>

#include "fcl/geometry/shape/tetrahedron.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/tetrahedron_intersect.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void benchmark_tetrahedron_with_tetrahedron() {
  const int tetrahedronCount = 1e3;
  const int test_n = 1e3;
  std::vector<Tetrahedron<S>> a_set;
  std::vector<Tetrahedron<S>> b_set;
  a_set.reserve(tetrahedronCount);
  b_set.reserve(tetrahedronCount);
  for (int i = 0; i < tetrahedronCount; i++) {
    a_set.emplace_back(test::generateRandomTetrahedron<S>());
    b_set.emplace_back(test::generateRandomTetrahedron<S>());
  }

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  std::vector<fcl::Transform3<S>> test_poses;
  test_poses.reserve(test_n);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> pose;
    fcl::test::generateRandomTransform(xyz_extent, pose);
    test_poses.emplace_back(pose);
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < tetrahedronCount; tet_index++) {
    const auto& a = a_set[tet_index];
    const auto& b = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      detail::tetrahedronTrahedronIntersect(a, Transform3<S>::Identity(), b,
                                            test_poses[test_idx]);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Tetrahedron_Tetrahedron_Time(by separating axis) in ms: "
            << ms_time << std::endl;

  start = std::chrono::high_resolution_clock::now();
  detail::GJKSolver<S> gjk_solver;
  for (int tet_index = 0; tet_index < tetrahedronCount; tet_index++) {
    const auto& a = a_set[tet_index];
    const auto& b = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      gjk_solver.shapeIntersect(a, Transform3<S>::Identity(), b,
                                test_poses[test_idx], nullptr);
    }
  }

  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "Tetrahedron_Tetrahedron_Time(by gjk) in ms: " << ms_time
            << std::endl;
}

template <typename S>
void benchmark_tetrahedron_with_box() {
  const int objectCount = 100;
  const int test_n = 1e4;
  std::vector<Tetrahedron<S>> a_set;
  std::vector<Box<S>> b_set;
  a_set.reserve(objectCount);
  b_set.reserve(objectCount);
  for (int i = 0; i < objectCount; i++) {
    a_set.emplace_back(test::generateRandomTetrahedron<S>());
    b_set.emplace_back(test::generateRandomBox<S>());
  }

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  std::vector<fcl::Transform3<S>> test_poses;
  test_poses.reserve(test_n);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> pose;
    fcl::test::generateRandomTransform(xyz_extent, pose);
    test_poses.emplace_back(std::move(pose));
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const Tetrahedron<S>& tetrahedron = a_set[tet_index];
    const auto& box = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      detail::boxTerahedronIntersect(box, Transform3<S>::Identity(),
                                     tetrahedron, test_poses[test_idx]);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Tetrahedron_Box_Time(by separating axis) in ms: " << ms_time
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  detail::GJKSolver<S> gjk_solver;
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const auto& a = a_set[tet_index];
    const auto& b = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      gjk_solver.shapeIntersect(a, Transform3<S>::Identity(), b,
                                test_poses[test_idx], nullptr);
    }
  }

  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "Tetrahedron_Box_Time(by gjk) in ms: " << ms_time << std::endl;
}

template <typename S>
void benchmark_tetrahedron_with_triangle() {
  const int objectCount = 100;
  const int test_n = 1e4;
  std::vector<Tetrahedron<S>> a_set;
  std::vector<TriangleP<S>> b_set;
  a_set.reserve(objectCount);
  b_set.reserve(objectCount);
  for (int i = 0; i < objectCount; i++) {
    a_set.emplace_back(test::generateRandomTetrahedron<S>());
    b_set.emplace_back(test::generateRandomTriangleP<S>());
  }

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  std::vector<fcl::Transform3<S>> test_poses;
  test_poses.reserve(test_n);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> pose;
    fcl::test::generateRandomTransform(xyz_extent, pose);
    test_poses.emplace_back(std::move(pose));
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const Tetrahedron<S>& tetrahedron = a_set[tet_index];
    const auto& tri = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      detail::triangleTerahedronIntersect(tri, Transform3<S>::Identity(),
                                          tetrahedron, test_poses[test_idx]);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Tetrahedron_Triangle_Time(by separating axis) in ms: "
            << ms_time << std::endl;

  detail::GJKSolver<S> gjk_solver;
  start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const auto& a = a_set[tet_index];
    const auto& b = b_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      gjk_solver.shapeIntersect(a, Transform3<S>::Identity(), b,
                                test_poses[test_idx], nullptr);
    }
  }

  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "Tetrahedron_Triangle_Time(by gjk) in ms: " << ms_time
            << std::endl;
}

std::shared_ptr<fcl::octomap::Pointcloud> generateRandomPointCloud(
    std::size_t numPoints) {
  auto cloud = std::make_shared<fcl::octomap::Pointcloud>();
  cloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    cloud->push_back(point[0], point[1], point[2] + 1);
  }
  return cloud;
}

template <typename BV>
void benchmark_tetrahedronBVHModel_with_heightmap() {
  int n_elem_x_axis = 10;
  int n_elem_y_axis = 10;
  using S = typename BV::S;
  Vector3<S> a{0, 0, 0};
  Vector3<S> b{1, 0, 0};
  Vector3<S> c{0, 0.5, 0};
  Vector3<S> d{0, 0, 1.5};
  const S x_offset = 1.0;
  const S y_offset = 2.0;

  // Start building
  BVHModel<BV> bvh;
  std::vector<Tetrahedron<S>> tetrahedrons;
  tetrahedrons.reserve(n_elem_x_axis * n_elem_y_axis);
  bvh.beginModel();

  // Add elements
  for (auto x_idx = 0; x_idx < n_elem_x_axis; x_idx++) {
    for (auto y_idx = 0; y_idx < n_elem_y_axis; y_idx++) {
      Vector3<S> displacement{x_offset * S(x_idx), y_offset * S(y_idx), 0};
      bvh.addTetrahedron(a + displacement, b + displacement, c + displacement,
                         d + displacement);
      tetrahedrons.emplace_back(a + displacement, b + displacement,
                                c + displacement, d + displacement);
    }
  }

  // Build the model
  bvh.endModel();
  const auto bvhPtr = std::make_shared<const BVHModel<BV>>(std::move(bvh));

  const auto points = generateRandomPointCloud(100000);
  auto heightMap =
      std::make_shared<fcl::heightmap::LayeredHeightMap<S>>(0.002, 512);
  heightMap->updateHeightsByPointCloud3D(*points);
  auto heightMap_geometry = fcl::HeightMapCollisionGeometry<S>(heightMap);
  heightMap_geometry.computeLocalAABB();

  detail::GJKSolver<S> narrowphase_solver;
  detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);

  // First checking
  const CollisionRequest<S> request{UINT_MAX};
  const int n_loop_per_pair = 50;
  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  const auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < n_loop_per_pair; i++) {
    Transform3<S> tf1;
    test::generateRandomTransform(extent, tf1);
    Transform3<S> tf2;
    test::generateRandomTransform(extent, tf2);

    CollisionResult<S> result;
    solver.HeightMapBVHIntersect(&heightMap_geometry, &bvh, tf1, tf2, request,
                                 result);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  const auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "TetrahedronBvh_HeightMap_Time in ms: " << ms_time << std::endl;
}

int main() {
  benchmark_tetrahedron_with_tetrahedron<double>();
  benchmark_tetrahedron_with_box<double>();
  benchmark_tetrahedron_with_triangle<double>();
  benchmark_tetrahedronBVHModel_with_heightmap<fcl::OBB<double>>();
  benchmark_tetrahedronBVHModel_with_heightmap<fcl::AABB<double>>();

  benchmark_tetrahedron_with_tetrahedron<float>();
  benchmark_tetrahedron_with_box<float>();
  benchmark_tetrahedron_with_triangle<float>();
  benchmark_tetrahedronBVHModel_with_heightmap<fcl::OBB<float>>();
  benchmark_tetrahedronBVHModel_with_heightmap<fcl::AABB<double>>();
  return 0;
}
