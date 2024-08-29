#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>

#include "fcl/cvx_collide/gjk.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_triangle.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void benchmark_box_with_triangle() {
  const int objectCount = 100;
  const int test_n = 1e4;
  std::vector<Box<S>> a_set;
  std::vector<TriangleP<S>> b_set;
  a_set.reserve(objectCount);
  b_set.reserve(objectCount);
  for (int i = 0; i < objectCount; i++) {
    a_set.emplace_back(test::generateRandomBox<S>());
    b_set.emplace_back(test::generateRandomTriangleP<S>());
  }

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  std::vector<fcl::Transform3<S>> test_poses(test_n);
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> pose;
    fcl::test::generateRandomTransform(xyz_extent, pose);
    test_poses.emplace_back(std::move(pose));
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const auto& tri = b_set[tet_index];
    const auto& box = a_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      detail::boxTriangleIntersect(box, Transform3<S>::Identity(), tri.a, tri.b,
                                   tri.c, test_poses[test_idx]);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Box_Triangle_Time(by separating axis) in ms: " << ms_time
            << std::endl;

  detail::GJKSolver<S> gjk_solver;
  start = std::chrono::high_resolution_clock::now();
  for (int tet_index = 0; tet_index < objectCount; tet_index++) {
    const auto& tri = b_set[tet_index];
    const auto& box = a_set[tet_index];
    for (int test_idx = 0; test_idx < test_n; test_idx++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&tri);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;
      minkowski_diff.toshape0 = test_poses[test_idx];
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();

      // Try with GJK/MPR
      detail::GJK2<S> gjk(1000, 1e-6);
      detail::GJKSimplex<S> simplex;
      gjk.Evaluate(minkowski_diff, simplex);
    }
  }

  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "Box_Triangle_Time(by gjk) in ms: " << ms_time << std::endl;
}

//==============================================================================
int main() {
  benchmark_box_with_triangle<double>();
  benchmark_box_with_triangle<float>();
  return 0;
}