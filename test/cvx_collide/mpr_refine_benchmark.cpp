//
// Created by mech-mind_gw on 3/29/2022.
//

#include <chrono>

#include "create_primitive_mesh.h"
#include "fcl/cvx_collide/epa.h"
#include "fcl/cvx_collide/mpr.h"
#include "retired_epa.h"
#include "retired_gjk.h"
#include "test_fcl_utility.h"

using namespace fcl;

template <typename S>
void minimumPenetrationDistanceBenchmark(
    const ShapeBase<S>* shape1, const ShapeBase<S>* shape2,
    std::function<bool(const fcl::Transform3<S>&, const fcl::Transform3<S>&,
                       fcl::Vector3<S>* direction, S* depth)>
        is_intersect_functor,
    const S disturb_angle_tangent = S(0.57),  // 30 degree
    const std::size_t test_n = 10000) {
  std::array<S, 6> extent{-0.3, -0.3, -0.3, 0.3, 0.3, 0.3};
  std::vector<fcl::Transform3<S>> shape2_pose_vector;
  std::vector<fcl::Vector3<S>> disturbed_direction_vector;
  std::vector<fcl::detail::GJK<S>> gjk_vector;
  std::vector<fcl::detail::MinkowskiDiff<S>> shape_vector;
  while (true) {
    fcl::Transform3<S> shape1_pose, shape2_pose;
    shape1_pose.setIdentity();
    test::generateRandomTransform(extent, shape2_pose);

    // Determine the intersection
    S gt_depth;
    Vector3<S> direction;
    bool is_intersect =
        is_intersect_functor(shape1_pose, shape2_pose, &direction, &gt_depth);
    if (!is_intersect) {
      continue;
    }

    // Direction must be valid
    if (direction.norm() <= 1e-3) {
      continue;
    }
    direction.normalize();

    // Make the shape
    detail::MinkowskiDiff<S> minkowski_diff;
    minkowski_diff.shapes[0] = detail::constructGJKGeometry(shape1);
    minkowski_diff.shapes[1] = detail::constructGJKGeometry(shape2);
    minkowski_diff.support_function = detail::computeSupport<S>;
    minkowski_diff.interior_function = detail::computeInterior<S>;
    minkowski_diff.toshape0 = shape1_pose.inverse() * shape2_pose;
    minkowski_diff.toshape1 =
        minkowski_diff.toshape0.inverse().rotation().matrix();
    shape_vector.emplace_back(std::move(minkowski_diff));

    // Run GJK/EPA
    detail::GJK<S> gjk(1000, 1e-6);
    Vector3<S> guess = direction;
    auto gjk_result = gjk.evaluate(shape_vector.back(), guess);
    if (gjk_result != detail::GJK<S>::Inside) {
      continue;
    }

    // Change the direction a
    Vector3<S> disturbed_direction;
    test::disturbDirectionByAngle(direction, disturb_angle_tangent,
                                  disturbed_direction);

    // Update the pose vector
    disturbed_direction_vector.emplace_back(std::move(disturbed_direction));
    shape2_pose_vector.emplace_back(std::move(shape2_pose));
    gjk_vector.emplace_back(std::move(gjk));
    if (shape2_pose_vector.size() >= test_n) {
      break;
    }
  }

  // Logging
  std::cout << "Data generation for " << test_n << ". Now test EPA"
            << std::endl;

  // Start EPA
  {
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < gjk_vector.size(); i++) {
      detail::GJK<S>& gjk_i = gjk_vector[i];
      detail::EPA<S> epa(256, 256, 1000, 1e-6);
      typename detail::EPA<S>::Status epa_status =
          epa.evaluate(gjk_i, disturbed_direction_vector[i]);
      if (epa_status != detail::EPA<S>::Failed) {
        // std::cout << "The epa depth " << epa.depth << std::endl;
      }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "EPA Time in ms " << ms_time << std::endl;
  }

  // Start EPA2
  {
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < gjk_vector.size(); i++) {
      detail::GJK<S>& gjk_i = gjk_vector[i];
      detail::EPA2<S> epa(128, 128, 1e-6);
      S depth;
      Vector3<S> p1, p2;
      // typename detail::EPA_Status epa_status =
      //     epa.Evaluate(gjk_i, gjk_i.shape, &depth, &p1, &p2);
      auto epa_status = fcl::detail::runEPA2_WithOldGJK(epa, gjk_i, gjk_i.shape,
                                                        &depth, &p1, &p2);

      if (epa_status != detail::EPA_Status::Failed) {
        // std::cout << "The epa depth " << depth << std::endl;
      }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "EPA2 Time in ms " << ms_time << std::endl;
  }

  // Start mpr
  {
    auto start = std::chrono::high_resolution_clock::now();
    detail::MPR<S> mpr(1000, 1e-6);
    for (std::size_t i = 0; i < gjk_vector.size(); i++) {
      typename detail::MPR<S>::IncrementalMinimumPenetrationData data;
      auto status = mpr.IncrementalMinimumPenetrationDistance(
          shape_vector[i], disturbed_direction_vector[i], data);
      assert(status == detail::MPR<S>::IncrementalPenetrationStatus::OK);
      FCL_UNUSED(status);  // disable warning in release config
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "MPR Time in ms " << ms_time << std::endl;
  }
}

template <typename S>
void spherePenetrationBenchmark() {
  fcl::Sphere<S> sphere1(S(0.15));
  fcl::Sphere<S> sphere2(S(0.15));
  auto intersect_functor = [](const fcl::Transform3<S>& pose1,
                              const fcl::Transform3<S>& pose2,
                              Vector3<S>* direction, S* depth) -> bool {
    fcl::Vector3<S> center_diff = pose2.translation() - pose1.translation();
    if (direction != nullptr) {
      *direction = center_diff.normalized();
    }
    if (depth != nullptr) {
      *depth = 0.3 - center_diff.norm();
    }
    return center_diff.norm() < 0.3;
  };
  minimumPenetrationDistanceBenchmark<S>(&sphere1, &sphere2, intersect_functor);
}

template <typename S>
void boxPenetrationBenchmark() {
  fcl::Box<S> box1(0.15, 0.15, 0.15);
  fcl::Box<S> box2(0.15, 0.15, 0.15);
  auto intersect_functor = [&box1, &box2](const fcl::Transform3<S>& pose1,
                                          const fcl::Transform3<S>& pose2,
                                          Vector3<S>* direction,
                                          S* depth) -> bool {
    CollisionPenetrationContactData<S> contacts;
    int return_code;
    Vector3<S> local_normal;
    S local_depth;
    detail::boxBox2(box1.side, pose1, box2.side, pose2, local_normal,
                    &local_depth, &return_code, 4, contacts);
    if (depth) *depth = local_depth;
    if (direction) *direction = local_normal;
    return return_code != 0;
  };

  minimumPenetrationDistanceBenchmark<S>(&box1, &box2, intersect_functor);
}

int main() {
  std::cout << "Benchmark sphere for double" << std::endl;
  spherePenetrationBenchmark<double>();
  std::cout << "Benchmark sphere for float" << std::endl;
  spherePenetrationBenchmark<float>();
  std::cout << "Benchmark box for double" << std::endl;
  boxPenetrationBenchmark<double>();
  std::cout << "Benchmark box for float" << std::endl;
  boxPenetrationBenchmark<float>();
}
