//
// Created by mech-mind_gw on 2021/12/29.
//

#include <array>
#include <chrono>
#include <limits>
#include <mutex>

#include "fcl/narrowphase/collision.h"
#include "fcl/math/detail/project.h"


class GeometryManager {
 public:
  void read() {
    std::shared_ptr<const fcl::CollisionGeometryf> my_copy;

    // Critical region
    {
      std::lock_guard<std::mutex> guard(mutex_);
      my_copy = geometry_;
    }

    // Do crazy staff on my_copy
  }

  void write() {
    // Prepare the geometry, can be very expensive
    std::shared_ptr<const fcl::CollisionGeometryf> my_copy;

    // Critical region
    {
      std::lock_guard<std::mutex> guard(mutex_);
      my_copy.swap(geometry_);
    }

    // Done
  }

 private:
  std::mutex mutex_;
  std::shared_ptr<const fcl::CollisionGeometryf> geometry_;
};

void triangle_intersection_test() {
  fcl::Boxf box(0.001, 0.001, 0.001);
  fcl::Spheref sphere(0.001);
  Eigen::Isometry3f pose;
  pose.setIdentity();
  pose.translation().x() = 0.1885;
  pose.translation().y() = -0.2755;
  pose.translation().z() = -0.0106;

  //
  Eigen::Vector3f p1, p2, p3;
  p1 = Eigen::Vector3f(0.185053, -0.275, 0.0334077);
  p2 = Eigen::Vector3f(0.0252288, -0.275, -0.186345);
  p3 = Eigen::Vector3f(0.187275, -0.275, -0.017001);

  // fcl::detail::GJKSolver_libccdf solver;
  fcl::detail::GJKSolver<float> solver;
  auto intersect_box = solver.shapeTriangleIntersect(box, pose, p1, p2, p3);
  auto intersect_sphere =
      solver.shapeTriangleIntersect(sphere, pose, p1, p2, p3);
  if (intersect_box && (!intersect_sphere))
    std::cout << "Box should be bounded by sphere" << std::endl;
}

void test_build_box_obbrss() {
  fcl::Boxf box(0.1, 0.3, 0.5);

  // Get a non-identity rotation
  float theta = 0.4;
  float c = std::cos(theta);
  float s = std::sin(theta);
  Eigen::Matrix3f rotation;
  rotation.setIdentity();
  rotation(1, 1) = c;
  rotation(1, 2) = s;
  rotation(2, 1) = -s;
  rotation(2, 2) = c;

  // Make the tf
  fcl::Transform3f box_tf;
  box_tf.setIdentity();
  box_tf.linear() = rotation;
  box_tf.translation() = Eigen::Vector3f(0.0, 0.1, 0.2);
  std::cout << box_tf.matrix() << std::endl;
  std::cout << rotation.transpose() * rotation << std::endl;

  // Make obbrssv
  int test_n = 1e6;
  auto start = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < test_n; i++) {
    fcl::OBBf bv;
    fcl::computeBV(box, box_tf, bv);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "Time is ms " << ms_time << std::endl;

  fcl::OBBRSSf bv_obbrss;
  fcl::computeBV(box, box_tf, bv_obbrss);
  std::cout << bv_obbrss.size() << std::endl;
}

template <typename T>
void test_project() {
  using namespace fcl::detail;
  fcl::Vector3<T> s1(0.000499996357,
                     -0.00349999964,
                     -0.00109999627);
  fcl::Vector3<T> s2(-0.000500003807,
                     0.0700000003,
                     0.00531465933);
  fcl::Vector3<T> s3(-0.000500003807,
                     -0.00349999964,
                     -0.00109999627);
  fcl::Vector3<T> s4(0.000499996357,
                     0.0700000003,
                     0.00531465933);
  //Vector3<T> v_min = parameterization[0] * s1 + parameterization[1] * s2 + parameterization[2] * s3 + parameterization[3] * s4;
  //std::cout << "v_min in test_project" << std::endl;
  //std::cout << v_min << std::endl;

  auto project_res = Project<T>::projectTetrahedraOrigin(s3, s2, s1, s4);
  // detail::Project<T>::ProjectResult project_res = detail::Project<T>::projectTriangleOrigin(s1, s2, s4);
  std::cout << "Old projection parameterization in test_project" << std::endl;
  for(auto i = 0; i < 4; i++) {
    std::cout << i << " " << project_res.parameterization[i] << std::endl;
  }
}

int main() {
  std::cout.precision(std::numeric_limits<float>::max_digits10);
  std::cout << "Project with float" << std::endl;
  test_project<float>();
  //std::cout << "Project with double" << std::endl;
  //test_project<double>();
  //triangle_intersection_test();
  // test_build_box_obbrss();
  return 0;
}