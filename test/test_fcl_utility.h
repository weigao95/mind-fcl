/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan */

#ifndef TEST_FCL_UTILITY_H
#define TEST_FCL_UTILITY_H

#include <gtest/gtest.h>

#include <array>
#include <fstream>
#include <iostream>

#include "fcl/common/unused.h"
#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/cylinder.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/math/constants.h"
#include "fcl/math/mesh_simplex.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/collision_object.h"
#include "fcl/narrowphase/collision_result.h"

#ifdef _WIN32
#define NOMINMAX  // required to avoid compilation errors with Visual Studio
                  // 2010
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace fcl {

namespace test {

class Timer {
 public:
  Timer();
  ~Timer();

  void start();                       ///< start timer
  void stop();                        ///< stop the timer
  double getElapsedTime();            ///< get elapsed time in milli-second
  double getElapsedTimeInSec();       ///< get elapsed time in second (same as
                                      ///< getElapsedTime)
  double getElapsedTimeInMilliSec();  ///< get elapsed time in milli-second
  double getElapsedTimeInMicroSec();  ///< get elapsed time in micro-second

 private:
  double startTimeInMicroSec;  ///< starting time in micro-second
  double endTimeInMicroSec;    ///< ending time in micro-second
  int stopped;                 ///< stop flag
#ifdef _WIN32
  LARGE_INTEGER frequency;  ///< ticks per second
  LARGE_INTEGER startCount;
  LARGE_INTEGER endCount;
#else
  timeval startCount;
  timeval endCount;
#endif
};

struct TStruct {
  std::vector<double> records;
  double overall_time;

  TStruct() { overall_time = 0; }

  void push_back(double t) {
    records.push_back(t);
    overall_time += t;
  }
};

/// @brief Load an obj mesh file
template <typename S>
void loadOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<MeshSimplex>& triangles);

template <typename S>
void saveOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<MeshSimplex>& triangles);

template <typename S>
S rand_interval(S rmin, S rmax);

template <typename S>
void eulerToMatrix(S a, S b, S c, Matrix3<S>& R);

template <typename S>
void disturbDirectionByAngle(const Vector3<S>& direction,
                             S disturb_angle_tangent, Vector3<S>& disturbed);

/// @brief Generate one random transform whose translation is constrained by
/// extents and rotation without constraints. The translation is (x, y, z), and
/// extents[0] <= x <= extents[3], extents[1] <= y <= extents[4], extents[2] <=
/// z <= extents[5]
template <typename S>
void generateRandomTransform(const std::array<S, 6>& extent,
                             Transform3<S>& transform);
template <typename S>
void generateRandomTransformVector(const std::array<S, 6>& extents,
                                   std::size_t n_samples,
                                   std::vector<Transform3<S>>& transform_vec);

/// @brief Generate n random transforms whose translations are constrained by
/// extents.
template <typename S>
void generateRandomTransforms(S extents[6],
                              aligned_vector<Transform3<S>>& transforms,
                              std::size_t n);

/// @brief Generate n random transforms whose translations are constrained by
/// extents. Also generate another transforms2 which have additional random
/// translation & rotation to the transforms generated.
template <typename S>
void generateRandomTransforms(S extents[6], S delta_trans[3], S delta_rot,
                              aligned_vector<Transform3<S>>& transforms,
                              aligned_vector<Transform3<S>>& transforms2,
                              std::size_t n);

/// @brief Generate n random tranforms and transform2 with addtional random
/// translation/rotation. The transforms and transform2 are used as initial and
/// goal configurations for the first mesh. The second mesh is in I. This is
/// used for continuous collision detection checking.
template <typename S>
void generateRandomTransforms_ccd(S extents[6],
                                  aligned_vector<Transform3<S>>& transforms,
                                  aligned_vector<Transform3<S>>& transforms2,
                                  S delta_trans[3], S delta_rot, std::size_t n,
                                  const std::vector<Vector3<S>>& vertices1,
                                  const std::vector<MeshSimplex>& triangles1,
                                  const std::vector<Vector3<S>>& vertices2,
                                  const std::vector<MeshSimplex>& triangles2);

/// @brief Generate environment with 3 * n objects: n boxes, n spheres and n
/// cylinders.
template <typename S>
void generateEnvironments(std::vector<CollisionObject<S>*>& env, S env_scale,
                          std::size_t n);

/// @brief Generate environment with 3 * n objects, but all in meshes.
template <typename S>
void generateEnvironmentsMesh(std::vector<CollisionObject<S>*>& env,
                              S env_scale, std::size_t n);

/// @brief Structure for minimum distance between two meshes and the
/// corresponding nearest point pair
template <typename S>
struct DistanceRes {
  S distance;
  Vector3<S> p1;
  Vector3<S> p2;
};

///
std::string getNodeTypeName(NODE_TYPE node_type);

/// @brief generate a bvh representing the given box
template <typename BV>
std::shared_ptr<const BVHModel<BV>> generateBoxBVHModel(
    const Box<typename BV::S>& box);

/// @brief make a octree2 with a plane
template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>> makePlane_xOy_AsOctree2(
    std::uint16_t bottom_half_shape, S scalar_resolution);
template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>> makeRandomPointsAsOctrees(
    S scalar_resolution, std::uint16_t bottom_half_shape,
    std::size_t n_random_points);

/// The directory of fcl/test
bool getTestModelDirectory(std::string& test_model_dir);

/// Build an octree from xyz point cloud file
std::shared_ptr<const octomap::Pointcloud> readPointCloudXYZ(
    const std::string& pcl_xyz_path);
std::shared_ptr<Octree2CollisionGeometry<float>> buildOctree2(
    const std::string& pcl_xyz_path, float voxel_resolution = 0.001,
    std::uint16_t bottom_half_shape = 1024);

/// Load a stl mesh from given path
void loadSTLMesh(const std::string& mesh_path,
                 std::vector<fcl::Vector3f>& points,
                 std::vector<fcl::MeshSimplex>& triangles);

/// Load a binvox file from given path
std::shared_ptr<fcl::octree2::Octree<float>> loadBinvoxAsOctree2(
    const std::string& filepath, double scale, std::uint16_t bottom_half_shape);

/// Printer
template <typename T>
std::string intersectStatusToString(
    typename fcl::detail::MPR<T>::IntersectStatus status) {
  using namespace fcl::detail;
  switch (status) {
    case MPR<T>::IntersectStatus::Intersect:
      return "Intersect";
    case MPR<T>::IntersectStatus::Separated:
      return "Separated";
    default:
      return "Failed";
  }
}

//============================================================================//
//                                                                            //
//                              Implementations                               //
//                                                                            //
//============================================================================//

//==============================================================================
template <typename S>
void loadOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<MeshSimplex>& triangles) {
  FILE* file = fopen(filename, "rb");
  if (!file) {
    std::cerr << "file not exist" << std::endl;
    return;
  }

  bool has_normal = false;
  bool has_texture = false;
  char line_buffer[2000];
  while (fgets(line_buffer, 2000, file)) {
    char* first_token = strtok(line_buffer, "\r\n\t ");
    if (!first_token || first_token[0] == '#' || first_token[0] == 0) continue;

    switch (first_token[0]) {
      case 'v': {
        if (first_token[1] == 'n') {
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          has_normal = true;
        } else if (first_token[1] == 't') {
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          has_texture = true;
        } else {
          S x = (S)atof(strtok(nullptr, "\t "));
          S y = (S)atof(strtok(nullptr, "\t "));
          S z = (S)atof(strtok(nullptr, "\t "));
          points.emplace_back(x, y, z);
        }
      } break;
      case 'f': {
        MeshSimplex tri;
        char* data[30];
        int n = 0;
        while ((data[n] = strtok(nullptr, "\t \r\n")) != nullptr) {
          if (strlen(data[n])) n++;
        }

        for (int t = 0; t < (n - 2); ++t) {
          if ((!has_texture) && (!has_normal)) {
            tri[0] = atoi(data[0]) - 1;
            tri[1] = atoi(data[1]) - 1;
            tri[2] = atoi(data[2]) - 1;
          } else {
            const char* v1;
            for (int i = 0; i < 3; i++) {
              // vertex ID
              if (i == 0)
                v1 = data[0];
              else
                v1 = data[t + i];

              tri[i] = atoi(v1) - 1;
            }
          }
          triangles.push_back(tri);
        }
      }
    }
  }
}

//==============================================================================
template <typename S>
void saveOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<MeshSimplex>& triangles) {
  std::ofstream os(filename);
  if (!os) {
    std::cerr << "file not exist" << std::endl;
    return;
  }

  for (std::size_t i = 0; i < points.size(); ++i) {
    os << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2]
       << std::endl;
  }

  for (std::size_t i = 0; i < triangles.size(); ++i) {
    os << "f " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << " "
       << triangles[i][2] + 1 << std::endl;
  }

  os.close();
}

//==============================================================================
template <typename S>
S rand_interval(S rmin, S rmax) {
  S t = rand() / ((S)RAND_MAX + 1);
  return (t * (rmax - rmin) + rmin);
}

//==============================================================================
template <typename S>
void eulerToMatrix(S a, S b, S c, Matrix3<S>& R) {
  auto c1 = std::cos(a);
  auto c2 = std::cos(b);
  auto c3 = std::cos(c);
  auto s1 = std::sin(a);
  auto s2 = std::sin(b);
  auto s3 = std::sin(c);

  R << c1 * c2, -c2 * s1, s2, c3 * s1 + c1 * s2 * s3, c1 * c3 - s1 * s2 * s3,
      -c2 * s3, s1 * s3 - c1 * c3 * s2, c3 * s1 * s2 + c1 * s3, c2 * c3;
}

//==============================================================================
template <typename S>
void disturbDirectionByAngle(const Vector3<S>& direction,
                             S disturb_angle_tangent, Vector3<S>& disturbed) {
  Vector3<S> random_disturb;
  while (true) {
    random_disturb.setRandom();
    if (random_disturb.squaredNorm() <= std::numeric_limits<S>::epsilon()) {
      continue;
    }

    // Compute the unit disturb perpendicular to the direction
    random_disturb.normalize();
    S dot_disturb_direction = random_disturb.dot(direction);
    random_disturb -= dot_disturb_direction * direction;
    // dot(random_disturb, direction) = 0 at this moment
    if (random_disturb.squaredNorm() <= std::numeric_limits<S>::epsilon()) {
      continue;
    }

    random_disturb.normalize();
    assert(std::abs(random_disturb.dot(direction)) <= 1e-4);
    disturbed = direction + disturb_angle_tangent * random_disturb;
    if (disturbed.squaredNorm() <= std::numeric_limits<S>::epsilon()) {
      continue;
    } else {
      disturbed.normalize();
      return;
    }
    // A direction found
  }
}

//==============================================================================
template <typename S>
void generateRandomTransform(const std::array<S, 6>& extents,
                             Transform3<S>& transform) {
  auto x = rand_interval(extents[0], extents[3]);
  auto y = rand_interval(extents[1], extents[4]);
  auto z = rand_interval(extents[2], extents[5]);

  const auto pi = constants<S>::pi();
  auto a = rand_interval((S)0, 2 * pi);
  auto b = rand_interval((S)0, 2 * pi);
  auto c = rand_interval((S)0, 2 * pi);

  Matrix3<S> R;
  eulerToMatrix(a, b, c, R);
  Vector3<S> T(x, y, z);
  transform.setIdentity();
  transform.linear() = R;
  transform.translation() = T;
}

//==============================================================================
template <typename S>
void generateRandomTransformVector(const std::array<S, 6>& extents,
                                   std::size_t n_samples,
                                   std::vector<Transform3<S>>& transform_vec) {
  transform_vec.resize(n_samples);
  for (std::size_t i = 0; i < n_samples; i++) {
    Transform3<S> sample_i;
    generateRandomTransform(extents, sample_i);
    transform_vec[i] = std::move(sample_i);
  }
}

//==============================================================================
template <typename S>
void generateRandomTransforms(S extents[6],
                              aligned_vector<Transform3<S>>& transforms,
                              std::size_t n) {
  transforms.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto x = rand_interval(extents[0], extents[3]);
    auto y = rand_interval(extents[1], extents[4]);
    auto z = rand_interval(extents[2], extents[5]);

    const auto pi = constants<S>::pi();
    auto a = rand_interval((S)0, 2 * pi);
    auto b = rand_interval((S)0, 2 * pi);
    auto c = rand_interval((S)0, 2 * pi);

    {
      Matrix3<S> R;
      eulerToMatrix(a, b, c, R);
      Vector3<S> T(x, y, z);
      transforms[i].setIdentity();
      transforms[i].linear() = R;
      transforms[i].translation() = T;
    }
  }
}

//==============================================================================
template <typename S>
void generateRandomTransforms(S extents[6], S delta_trans[3], S delta_rot,
                              aligned_vector<Transform3<S>>& transforms,
                              aligned_vector<Transform3<S>>& transforms2,
                              std::size_t n) {
  transforms.resize(n);
  transforms2.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto x = rand_interval(extents[0], extents[3]);
    auto y = rand_interval(extents[1], extents[4]);
    auto z = rand_interval(extents[2], extents[5]);

    const auto pi = constants<S>::pi();
    auto a = rand_interval((S)0, 2 * pi);
    auto b = rand_interval((S)0, 2 * pi);
    auto c = rand_interval((S)0, 2 * pi);

    {
      Matrix3<S> R;
      eulerToMatrix(a, b, c, R);
      Vector3<S> T(x, y, z);
      transforms[i].setIdentity();
      transforms[i].linear() = R;
      transforms[i].translation() = T;
    }

    auto deltax = rand_interval(-delta_trans[0], delta_trans[0]);
    auto deltay = rand_interval(-delta_trans[1], delta_trans[1]);
    auto deltaz = rand_interval(-delta_trans[2], delta_trans[2]);

    auto deltaa = rand_interval(-delta_rot, delta_rot);
    auto deltab = rand_interval(-delta_rot, delta_rot);
    auto deltac = rand_interval(-delta_rot, delta_rot);

    {
      Matrix3<S> R;
      eulerToMatrix(a + deltaa, b + deltab, c + deltac, R);
      Vector3<S> T(x + deltax, y + deltay, z + deltaz);
      transforms2[i].setIdentity();
      transforms2[i].linear() = R;
      transforms2[i].translation() = T;
    }
  }
}

//==============================================================================
template <typename S>
void generateRandomTransforms_ccd(S extents[6],
                                  aligned_vector<Transform3<S>>& transforms,
                                  aligned_vector<Transform3<S>>& transforms2,
                                  S delta_trans[3], S delta_rot, std::size_t n,
                                  const std::vector<Vector3<S>>& vertices1,
                                  const std::vector<MeshSimplex>& triangles1,
                                  const std::vector<Vector3<S>>& vertices2,
                                  const std::vector<MeshSimplex>& triangles2) {
  FCL_UNUSED(vertices1);
  FCL_UNUSED(triangles1);
  FCL_UNUSED(vertices2);
  FCL_UNUSED(triangles2);

  transforms.resize(n);
  transforms2.resize(n);

  for (std::size_t i = 0; i < n;) {
    auto x = rand_interval(extents[0], extents[3]);
    auto y = rand_interval(extents[1], extents[4]);
    auto z = rand_interval(extents[2], extents[5]);

    const auto pi = constants<S>::pi();
    auto a = rand_interval(0, 2 * pi);
    auto b = rand_interval(0, 2 * pi);
    auto c = rand_interval(0, 2 * pi);

    Matrix3<S> R;
    eulerToMatrix(a, b, c, R);
    Vector3<S> T(x, y, z);
    Transform3<S> tf(Transform3<S>::Identity());
    tf.linear() = R;
    tf.translation() = T;

    std::vector<std::pair<int, int>> results;
    {
      transforms[i] = tf;

      auto deltax = rand_interval(-delta_trans[0], delta_trans[0]);
      auto deltay = rand_interval(-delta_trans[1], delta_trans[1]);
      auto deltaz = rand_interval(-delta_trans[2], delta_trans[2]);

      auto deltaa = rand_interval(-delta_rot, delta_rot);
      auto deltab = rand_interval(-delta_rot, delta_rot);
      auto deltac = rand_interval(-delta_rot, delta_rot);

      Matrix3<S> R2;
      eulerToMatrix(a + deltaa, b + deltab, c + deltac, R2);
      Vector3<S> T2(x + deltax, y + deltay, z + deltaz);
      transforms2[i].linear() = R2;
      transforms2[i].translation() = T2;
      ++i;
    }
  }
}

//==============================================================================
template <typename S>
void generateEnvironments(std::vector<CollisionObject<S>*>& env, S env_scale,
                          std::size_t n) {
  S extents[] = {-env_scale, env_scale,  -env_scale,
                 env_scale,  -env_scale, env_scale};
  aligned_vector<Transform3<S>> transforms;

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Box<S>* box = new Box<S>(5, 10, 20);
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(box), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Sphere<S>* sphere = new Sphere<S>(30);
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(sphere), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Cylinder<S>* cylinder = new Cylinder<S>(10, 40);
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(cylinder), transforms[i]));
  }
}

//==============================================================================
template <typename S>
void generateEnvironmentsMesh(std::vector<CollisionObject<S>*>& env,
                              S env_scale, std::size_t n) {
  S extents[] = {-env_scale, env_scale,  -env_scale,
                 env_scale,  -env_scale, env_scale};
  aligned_vector<Transform3<S>> transforms;

  generateRandomTransforms(extents, transforms, n);
  Box<S> box(5, 10, 20);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, box, Transform3<S>::Identity());
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  Sphere<S> sphere(30);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, sphere, Transform3<S>::Identity(), 16, 16);
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  Cylinder<S> cylinder(10, 40);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, cylinder, Transform3<S>::Identity(), 16, 16);
    env.push_back(new CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }
}

//==============================================================================
template <typename BV>
std::shared_ptr<const BVHModel<BV>> generateBoxBVHModel(
    const Box<typename BV::S>& box) {
  using S = typename BV::S;
  // Box bvh
  std::shared_ptr<BVHModel<BV>> model(new BVHModel<BV>());
  auto a = box.side[0];
  auto b = box.side[1];
  auto c = box.side[2];

  std::vector<Vector3<S>> points(8);
  std::vector<MeshSimplex> tri_indices(12);
  points[0] << 0.5 * a, -0.5 * b, 0.5 * c;
  points[1] << 0.5 * a, 0.5 * b, 0.5 * c;
  points[2] << -0.5 * a, 0.5 * b, 0.5 * c;
  points[3] << -0.5 * a, -0.5 * b, 0.5 * c;
  points[4] << 0.5 * a, -0.5 * b, -0.5 * c;
  points[5] << 0.5 * a, 0.5 * b, -0.5 * c;
  points[6] << -0.5 * a, 0.5 * b, -0.5 * c;
  points[7] << -0.5 * a, -0.5 * b, -0.5 * c;

  tri_indices[0].set(0, 4, 1);
  tri_indices[1].set(1, 4, 5);
  tri_indices[2].set(2, 6, 3);
  tri_indices[3].set(3, 6, 7);
  tri_indices[4].set(3, 0, 2);
  tri_indices[5].set(2, 0, 1);
  tri_indices[6].set(6, 5, 7);
  tri_indices[7].set(7, 5, 4);
  tri_indices[8].set(1, 5, 2);
  tri_indices[9].set(2, 5, 6);
  tri_indices[10].set(3, 7, 0);
  tri_indices[11].set(0, 7, 4);

  int result;
  result = model->beginModel();
  EXPECT_EQ(result, BVH_OK);

  for (std::size_t i = 0; i < tri_indices.size(); ++i) {
    result =
        model->addTriangle(points[tri_indices[i][0]], points[tri_indices[i][1]],
                           points[tri_indices[i][2]]);
    EXPECT_EQ(result, BVH_OK);
  }

  result = model->endModel();
  EXPECT_EQ(result, BVH_OK);
  model->computeLocalAABB();

  EXPECT_EQ(model->num_vertices(), 12 * 3);
  EXPECT_EQ(model->num_simplex(), 12);
  EXPECT_EQ(model->build_state(), BVH_BUILD_STATE_PROCESSED);
  return model;
}

//==============================================================================
template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>> makePlane_xOy_AsOctree2(
    std::uint16_t bottom_half_shape, S scalar_resolution) {
  std::vector<Vector3<S>> points;
  for (int i_x = 0; i_x < bottom_half_shape; i_x++) {
    for (int i_y = 0; i_y < bottom_half_shape; i_y++) {
      Vector3<S> point;
      point.x() = (S(i_x) + 0.5) * scalar_resolution;
      point.y() = (S(i_y) + 0.5) * scalar_resolution;
      point.z() = 0.1 * scalar_resolution;
      points.push_back(point);
      point.z() = -0.1 * scalar_resolution;
      points.push_back(point);
    }
  }

  Vector3<S> resolution(scalar_resolution, scalar_resolution,
                        scalar_resolution);
  auto tree =
      std::make_shared<octree2::Octree<S>>(resolution, bottom_half_shape);

  auto point_fn = [&points](int index, S& x, S& y, S& z) -> void {
    const auto point = points[index];
    x = point.x();
    y = point.y();
    z = point.z();
  };
  tree->rebuildTree(point_fn, points.size());
  auto tree_geom = std::make_shared<Octree2CollisionGeometry<S>>(tree);
  tree_geom->computeLocalAABB();
  return tree_geom;
}

//==============================================================================
template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>> makeRandomPointsAsOctrees(
    S scalar_resolution, std::uint16_t bottom_half_shape,
    std::size_t n_random_points) {
  using namespace octree2;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  Vector3<S> resolution(scalar_resolution, scalar_resolution,
                        scalar_resolution);
  std::vector<Vector3<S>> points_inserted;
  for (std::size_t i = 0; i < n_random_points; i++) {
    Vector3<S> point_i;
    point_i.setRandom();
    point_i *= (0.99 * bottom_half_size);
    points_inserted.push_back(point_i);
  }

  // Insert into tree
  auto point_fn = [&points_inserted](int index, S& x, S& y, S& z) -> void {
    const auto point = points_inserted[index];
    x = point.x();
    y = point.y();
    z = point.z();
  };
  auto tree = std::make_shared<Octree<S>>(resolution, bottom_half_shape);
  tree->rebuildTree(point_fn, points_inserted.size());

  auto octree_geom = std::make_shared<Octree2CollisionGeometry<S>>(tree);
  octree_geom->computeLocalAABB();
  return octree_geom;
}

//==============================================================================
template <typename S>
Tetrahedron<S> generateRandomTetrahedron(S min_size = 0.02, S max_size = 0.2) {
  const Vector3<S> v0{0., 0., 0.};
  const Vector3<S> v1{fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size)};
  const Vector3<S> v2{fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size)};
  const Vector3<S> v3{fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size)};
  Tetrahedron<S> tet(v0, v1, v2, v3);
  tet.computeLocalAABB();
  return tet;
}

template <typename S>
Box<S> generateRandomBox(S min_size = 0.02, S max_size = 0.2) {
  const S x = fcl::test::rand_interval<S>(min_size, max_size);
  const S y = fcl::test::rand_interval<S>(min_size, max_size);
  const S z = fcl::test::rand_interval<S>(min_size, max_size);
  Box<S> box{x, y, z};
  box.computeLocalAABB();
  return box;
}

template <typename S>
const TriangleP<S> generateRandomTriangleP(S min_size = 0.02,
                                           S max_size = 0.2) {
  const Vector3<S> v0{0., 0., 0.};
  const Vector3<S> v1{fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size)};
  const Vector3<S> v2{fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size),
                      fcl::test::rand_interval<S>(min_size, max_size)};
  TriangleP<S> triangle{v0, v1, v2};
  triangle.computeLocalAABB();
  return triangle;
}

}  // namespace test
}  // namespace fcl

#endif
