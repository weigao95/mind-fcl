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

#ifndef FCL_MATH_TRIANGLE_H
#define FCL_MATH_TRIANGLE_H

#include <array>
#include <cstddef>

#include "fcl/common/types.h"

namespace fcl {

/// @brief Triangle or tetrahedron with 3 or 4 indices for points
class MeshSimplex {
  /// @brief indices for each vertex of triangle
  std::array<int, 4> vertex_index{-1, -1, -1, -1};
  static constexpr int kInvalidVertexIndex = -1;

 public:
  /// @brief Default constructor
  MeshSimplex() = default;

  /// @brief Create a triangle with given vertex indices
  MeshSimplex(int point_idx_1, int point_idx_2, int point_idx_3) {
    set(point_idx_1, point_idx_2, point_idx_3);
  }

  /// @brief Create a tetrahedron with given vertex indices
  MeshSimplex(int point_idx_1, int point_idx_2, int point_idx_3,
              int point_idx_4) {
    set(point_idx_1, point_idx_2, point_idx_3, point_idx_4);
  }

  /// @brief Set the vertex indices of the simplex
  inline void set(int point_idx_1, int point_idx_2, int point_idx_3) {
    vertex_index[0] = point_idx_1;
    vertex_index[1] = point_idx_2;
    vertex_index[2] = point_idx_3;
    vertex_index[3] = kInvalidVertexIndex;
  }

  /// @brief Set the tetrahedron indices of the simplex
  inline void set(int point_idx_1, int point_idx_2, int point_idx_3,
                  int point_idx_4) {
    vertex_index[0] = point_idx_1;
    vertex_index[1] = point_idx_2;
    vertex_index[2] = point_idx_3;
    vertex_index[3] = point_idx_4;
  }

  /// @access the simplex index
  inline int operator[](int i) const { return vertex_index[i]; }
  inline int& operator[](int i) { return vertex_index[i]; }
  inline bool is_triangle() const { return vertex_index[3] < 0; }
  inline bool is_tetrahedron() const { return vertex_index[3] >= 0; }
};

template <typename S>
class Simplex {
 public:
  explicit Simplex() = default;
  Simplex(const Vector3<S>& a, const Vector3<S>& b, const Vector3<S>& c)
      : is_triangle_(true) {
    simplex_points_[0] = a;
    simplex_points_[1] = b;
    simplex_points_[2] = c;
  };
  Simplex(const Vector3<S>& a, const Vector3<S>& b, const Vector3<S>& c,
          const Vector3<S>& d)
      : is_triangle_(false) {
    simplex_points_[0] = a;
    simplex_points_[1] = b;
    simplex_points_[2] = c;
    simplex_points_[3] = d;
  };
  ~Simplex() = default;

  /// Scale the simplex wrt the center point
  void scaleWrtCenter(S scale) {
    // First compute the center
    Vector3<S> center;
    {
      const auto& p = simplex_points_;
      center = is_triangle_ ? ((p[0] + p[1] + p[2]).eval() / 3)
                            : ((p[0] + p[1] + p[2] + p[3]).eval() / 4);
    }

    // Then scale it
    for(auto i = 0; i < get_num_points(); i++) {
      Vector3<S> point_i = simplex_points_[i];
      Vector3<S> scaled_point_i = (point_i - center) * scale + center;
      simplex_points_[i] = scaled_point_i;
    }
  }

  /// @access the simplex index
  inline int get_num_points() const { return is_triangle_ ? 3 : 4; }
  inline bool is_triangle() const { return is_triangle_; }
  inline bool is_tetrahedron() const { return !is_triangle_; }
  inline Vector3<S>& operator[](int i) { return simplex_points_[i]; }
  inline const Vector3<S>& operator[](int i) const {
    return simplex_points_[i];
  }

 private:
  std::array<Vector3<S>, 4> simplex_points_{};
  bool is_triangle_{true};
};

}  // namespace fcl

#endif
