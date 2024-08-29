//
// Created by mech-mind_gw on 1/10/2023.
//

#include <chrono>

#include "create_primitive_mesh.h"
#include "retired_gjk.h"
#include "fcl/cvx_collide/gjk.h"
#include "fcl/cvx_collide/mpr.h"
#include "fcl/geometry/shape/shape_gjk_interface.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void triangleBoxCollisionBenchmark() {
  constexpr int test_n = 4e6;
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
  fcl::Box<S> box(0.4, 0.4, 0.4);

  // Test with gjk
  {
    detail::GJK<S> gjk(1000, 1e-6);
    auto start = std::chrono::high_resolution_clock::now();
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
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK Time in ms " << ms_time << std::endl;
  }

  // Test with gjk2
  {
    detail::GJK2<S> gjk2(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&box);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&triangles_vector[i]);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;
      minkowski_diff.toshape0.setIdentity();
      minkowski_diff.toshape1.setIdentity();
      gjk2.Evaluate(minkowski_diff, simplex, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK2 Time in ms " << ms_time << std::endl;
  }
}

template <typename S>
void capsuleCapsuleCollisionBenchmark() {
  fcl::Capsule<S> capsule_1(0.5, 1.3);
  fcl::Capsule<S> capsule_2(0.5, 1.3);
  int test_n = 1e5;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};

  // Make the pose
  std::vector<fcl::Transform3<S>> pose1_vector, pose2_vector;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);
    pose1_vector.emplace_back(obj1_pose);
    pose2_vector.emplace_back(obj2_pose);
  }

  // Test with gjk
  {
    detail::GJK<S> gjk(1000, 1e-6);
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK Time in ms " << ms_time << std::endl;
  }

  // Test with gjk2
  {
    detail::GJK2<S> gjk2(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk2.Evaluate(minkowski_diff, simplex, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK2 Time in ms " << ms_time << std::endl;
  }

  // Test with primitive
  {
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      // Setup the pose and run it
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      S capsule_distance;
      fcl::Vector3<S> point_1, point_2;
      detail::capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2,
                                        obj2_pose, &capsule_distance, &point_1,
                                        &point_2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "Primitive Time in ms " << ms_time << std::endl;
  }
}

template <typename S>
void capsuleSeparationBenchmark() {
  fcl::Capsule<S> capsule_1(0.5, 1.3);
  fcl::Capsule<S> capsule_2(0.5, 1.3);
  std::size_t test_n = 1e5;
  std::array<S, 6> extent{-2, -2, -2, 2, 2, 2};

  // Make the pose
  std::vector<fcl::Transform3<S>> pose1_vector, pose2_vector;
  while (pose1_vector.size() < test_n) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform(extent, obj1_pose);
    test::generateRandomTransform(extent, obj2_pose);

    // Do not emplace if not separated
    S capsule_distance;
    fcl::Vector3<S> point_1, point_2;
    detail::capsuleCapsuleDistance<S>(capsule_1, obj1_pose, capsule_2,
                                      obj2_pose, &capsule_distance, &point_1,
                                      &point_2);
    if (capsule_distance >= 0) {
      continue;
    } else {
      pose1_vector.emplace_back(obj1_pose);
      pose2_vector.emplace_back(obj2_pose);
    }
  }

  // Test with gjk
  {
    detail::GJK<S> gjk(1000, 1e-6);
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK Time in ms " << ms_time << std::endl;
  }

  // Test with gjk2
  {
    detail::GJK2<S> gjk2(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    typename detail::GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(&capsule_1);
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(&capsule_2);
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk2.Evaluate(minkowski_diff, simplex, Vector3<S>::UnitX(),
                    &min_distance_output);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK2 Time in ms " << ms_time << std::endl;
  }
}

template <typename S>
void convexSphereSeparationBenchmark() {
  const S sphere_radius = 1.0;
  const S expected_distance = 0.2;
  const S displacement = sphere_radius * 2 + expected_distance;
  auto sphere1 = std::make_shared<fcl::Convex<S>>(
      fcl::test::makeSphereConvex<S>(sphere_radius));
  auto sphere2 = std::make_shared<fcl::Convex<S>>(
      fcl::test::makeSphereConvex<S>(sphere_radius));

  // Construct the objects
  int test_n = 1e5;
  std::vector<fcl::Transform3<S>> pose1_vector, pose2_vector;
  std::array<S, 6> extent{-3, -3, -3, 3, 3, 3};
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> sphere1_pose, pose_tmp, sphere2_pose;
    test::generateRandomTransform(extent, pose_tmp);
    fcl::Vector3<S> direction_in_world =
        pose_tmp.matrix().template block<3, 1>(0, 0);
    test::generateRandomTransform(extent, sphere1_pose);
    test::generateRandomTransform(extent, sphere2_pose);
    sphere2_pose.translation() =
        sphere1_pose.translation() + direction_in_world * displacement;

    pose1_vector.emplace_back(sphere1_pose);
    pose2_vector.emplace_back(sphere2_pose);
  }

  // Test with gjk
  {
    detail::GJK<S> gjk(1000, 1e-6);
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(sphere1.get());
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(sphere2.get());
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK Time in ms " << ms_time << std::endl;
  }

  // Test with gjk2
  {
    detail::GJK2<S> gjk2(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    typename detail::GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(sphere1.get());
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(sphere2.get());
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk2.Evaluate(minkowski_diff, simplex, Vector3<S>::UnitX(),
                    &min_distance_output);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK2 Time in ms " << ms_time << std::endl;
  }
}

template <typename S>
void convexSphereBenchmark() {
  constexpr S radius = 0.1;
  constexpr int test_n = 1e6;
  auto sphere_1 =
      std::make_shared<fcl::Convex<S>>(fcl::test::makeSphereConvex<S>(radius));
  auto sphere_2 =
      std::make_shared<fcl::Convex<S>>(fcl::test::makeSphereConvex<S>(radius));

  // The extent used to generate the pose
  S extent = radius * 2;
  std::array<S, 6> xyz_extent{-extent, -extent, -extent,
                              extent,  extent,  extent};
  std::vector<fcl::Transform3<S>> pose1_vector, pose2_vector;
  for (auto test_idx = 0; test_idx < test_n; test_idx++) {
    fcl::Transform3<S> obj1_pose, obj2_pose;
    test::generateRandomTransform<S>(xyz_extent, obj1_pose);
    test::generateRandomTransform<S>(xyz_extent, obj2_pose);
    pose1_vector.emplace_back(obj1_pose);
    pose2_vector.emplace_back(obj2_pose);
  }

  // Test with gjk
  {
    detail::GJK<S> gjk(1000, 1e-6);
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(sphere_1.get());
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(sphere_2.get());
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk.evaluate(minkowski_diff, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK Time in ms " << ms_time << std::endl;
  }

  // Test with gjk2
  {
    detail::GJK2<S> gjk2(1000, 1e-6);
    detail::GJKSimplex<S> simplex;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < test_n; i++) {
      detail::MinkowskiDiff<S> minkowski_diff;
      minkowski_diff.shapes[0] = detail::constructGJKGeometry(sphere_1.get());
      minkowski_diff.shapes[1] = detail::constructGJKGeometry(sphere_2.get());
      minkowski_diff.support_function = detail::computeSupport<S>;
      minkowski_diff.interior_function = detail::computeInterior<S>;

      // Setup the pose
      const fcl::Transform3<S>& obj1_pose = pose1_vector[i];
      const fcl::Transform3<S>& obj2_pose = pose2_vector[i];
      minkowski_diff.toshape0 = obj1_pose.inverse() * obj2_pose;
      minkowski_diff.toshape1 =
          minkowski_diff.toshape0.inverse().rotation().matrix();
      gjk2.Evaluate(minkowski_diff, simplex, Vector3<S>::UnitX());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "GJK2 Time in ms " << ms_time << std::endl;
  }
}

void runTriangleBoxBenchmark() {
  std::cout << "TriangleBoxBenchmark with float" << std::endl;
  triangleBoxCollisionBenchmark<float>();
  std::cout << "TriangleBoxBenchmark with double" << std::endl;
  triangleBoxCollisionBenchmark<double>();
}

void runCapsuleCapsuleBenchmark() {
  std::cout << "CapsuleCapsuleBenchmark with float" << std::endl;
  capsuleCapsuleCollisionBenchmark<float>();
  std::cout << "CapsuleCapsuleBenchmark with double" << std::endl;
  capsuleCapsuleCollisionBenchmark<double>();
}

void runSphereConvexBenchmark() {
  std::cout << "SphereConvex with float" << std::endl;
  convexSphereBenchmark<float>();
  std::cout << "SphereConvex with double" << std::endl;
  convexSphereBenchmark<double>();
}

void runCapsuleSeparationBenchmark() {
  std::cout << "CapsuleSeparationBenchmark with float" << std::endl;
  capsuleSeparationBenchmark<float>();
  std::cout << "CapsuleSeparationBenchmark with double" << std::endl;
  capsuleSeparationBenchmark<double>();
}

void runSphereConvexSeparationDistanceBenchmark() {
  std::cout << "SphereConvexSeparationDistance with float" << std::endl;
  convexSphereSeparationBenchmark<float>();
  std::cout << "SphereConvexSeparationDistance with double" << std::endl;
  convexSphereSeparationBenchmark<double>();
}

int main() {
  runTriangleBoxBenchmark();
  runCapsuleCapsuleBenchmark();
  runSphereConvexBenchmark();
  runCapsuleSeparationBenchmark();
  runSphereConvexSeparationDistanceBenchmark();
}
