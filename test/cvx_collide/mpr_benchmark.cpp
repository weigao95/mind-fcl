//
// Created by mech-mind_gw on 2/8/2022.
//

#include <chrono>

#include "retired_gjk.h"
#include "fcl/narrowphase/detail/gjk_solver_cvx.h"
#include "fcl/cvx_collide/mpr.h"

using namespace fcl;

template <typename S>
void triangleBoxBenchmark() {
  constexpr int test_n = 1e6;
  std::vector<fcl::TriangleP<S>> triangles_vector;
  triangles_vector.reserve(test_n);
  for (auto i = 0; i < test_n; i++) {
    // A random triangle
    Vector3<S> a, b, c;
    a.setRandom();
    b.setRandom();
    c.setRandom();
    fcl::TriangleP<S> triangle(a, b, c);
    triangles_vector.emplace_back(std::move(triangle));
  }

  // Make the minkowski difference
  fcl::Box<S> box(0.1, 0.1, 0.1);

  // Test the mpr
  detail::MPR<S> mpr(1000, 1e-6);
  typename detail::MPR<S>::IntersectData intersect_data;
  auto start = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < test_n; i++) {
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangles_vector[i]);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0.setIdentity();
    minkowski_diff.toshape1.setIdentity();
    mpr.Intersect(minkowski_diff, &intersect_data);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "MPR Time in ms " << ms_time << std::endl;

  // Test with gjk
  detail::GJK<S> gjk(1000, 1e-6);
  start = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < test_n; i++) {
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangles_vector[i]);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0.setIdentity();
    minkowski_diff.toshape1.setIdentity();
    gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
  }
  end = std::chrono::high_resolution_clock::now();
  ms_time = std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
  std::cout << "GJK Time in ms " << ms_time << std::endl;
}

int main() {
  std::cout << "Test with float" << std::endl;
  triangleBoxBenchmark<float>();
  std::cout << "Test with double" << std::endl;
  triangleBoxBenchmark<double>();
}
