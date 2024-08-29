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

#ifndef FCL_NARROWPHASE_GJKSOLVERINDEP_H
#define FCL_NARROWPHASE_GJKSOLVERINDEP_H

#include <iostream>

#include "fcl/common/types.h"
#include "fcl/narrowphase/contact_point.h"
#include "fcl/narrowphase/detail/collision_penetration_mode.h"

namespace fcl {

namespace detail {

/// @brief collision and distance solver based on GJK algorithm implemented in
/// fcl (rewritten the code from the GJK in bullet)
template <typename S_>
struct GJKSolver {
  using S = S_;

  /// @brief intersection checking between two shapes
  template <typename Shape1, typename Shape2>
  bool shapeIntersect(
      const Shape1& s1, const Transform3<S>& tf1, const Shape2& s2,
      const Transform3<S>& tf2,
      CollisionPenetrationContactData<S>* penetration_contacts = nullptr) const;

  /// @brief intersection checking between one shape and a triangle
  template <typename Shape>
  bool shapeTriangleIntersect(
      const Shape& s, const Transform3<S>& tf, const Vector3<S>& P1,
      const Vector3<S>& P2, const Vector3<S>& P3,
      CollisionPenetrationContactData<S>* penetration_contacts = nullptr) const;

  //// @brief intersection checking between one shape and a triangle with
  /// transformation
  template <typename Shape>
  bool shapeTriangleIntersect(
      const Shape& s, const Transform3<S>& tf1, const Vector3<S>& P1,
      const Vector3<S>& P2, const Vector3<S>& P3, const Transform3<S>& tf2,
      CollisionPenetrationContactData<S>* penetration_contacts = nullptr) const;

  //// @brief intersection checking between one shape and a tetrahedron with
  /// transformation
  template <typename Shape>
  bool shapeTetrahedronIntersect(
      const Shape& s, const Transform3<S>& tf1, const Vector3<S>& P1,
      const Vector3<S>& P2, const Vector3<S>& P3, const Vector3<S>& P4,
      const Transform3<S>& tf2,
      CollisionPenetrationContactData<S>* penetration_contacts = nullptr) const;

  /// @brief distance computation between two shapes
  template <typename Shape1, typename Shape2>
  bool shapeDistance(const Shape1& s1, const Transform3<S>& tf1,
                     const Shape2& s2, const Transform3<S>& tf2,
                     S* distance = nullptr, Vector3<S>* p1 = nullptr,
                     Vector3<S>* p2 = nullptr) const;

  /// @brief distance computation between two shapes
  template <typename Shape1, typename Shape2>
  bool shapeSignedDistance(const Shape1& s1, const Transform3<S>& tf1,
                           const Shape2& s2, const Transform3<S>& tf2,
                           S* distance = nullptr, Vector3<S>* p1 = nullptr,
                           Vector3<S>* p2 = nullptr) const;

  /// @brief distance computation between one shape and a triangle
  template <typename Shape>
  bool shapeTriangleDistance(const Shape& s, const Transform3<S>& tf,
                             const Vector3<S>& P1, const Vector3<S>& P2,
                             const Vector3<S>& P3, S* distance = nullptr,
                             Vector3<S>* p1 = nullptr,
                             Vector3<S>* p2 = nullptr) const;

  /// @brief distance computation between one shape and a triangle with
  /// transformation
  template <typename Shape>
  bool shapeTriangleDistance(const Shape& s, const Transform3<S>& tf1,
                             const Vector3<S>& P1, const Vector3<S>& P2,
                             const Vector3<S>& P3, const Transform3<S>& tf2,
                             S* distance = nullptr, Vector3<S>* p1 = nullptr,
                             Vector3<S>* p2 = nullptr) const;

  /// @brief default setting for GJK algorithm
  GJKSolver();

  /// @brief maximum number of simplex face used in EPA algorithm
  unsigned int epa_max_face_num;

  /// @brief maximum number of simplex vertex used in EPA algorithm
  unsigned int epa_max_vertex_num;

  /// @brief maximum number of iterations used for EPA iterations
  unsigned int epa_max_iterations;

  /// @brief the threshold used in EPA to stop iteration
  S epa_tolerance;

  /// @brief the threshold used in GJK to stop iteration
  S gjk_tolerance;

  /// @brief maximum number of iterations used for GJK iterations
  S gjk_max_iterations;

  /// @brief Whether a direction guess can be provided for GJK
  ///        If two objects are disjoint, then a direction that the GJK
  ///        algorithm can determine the MinkowskiDiff are separated
  ///        would be beneficial to
  ///        1. speedup the determination of separation
  ///        2. If the provided direction are close to minimum separation
  ///           direction, and separation distance is requested. Then this
  ///           direction also speed up separation distance computation.
  ///        From the discussion above, the direction guess should be selected
  ///        as the min separation direction in SHAPE_1 Frame.
  ///
  ///        Something to Note about this direction hint:
  ///        1. EPA and MPR do NOT take this hint as input
  ///        2. If the GJKSolver is used as a sub-step of traversal algorithm,
  ///           such as octree collision detection. Then it would be very hard
  ///           to guess the separation direction, as collision detection are
  ///           performed with boxes in the octree.
  ///
  ///        As mentioned above, this direction hint is useful only for rather
  ///        limited case. Thus, we plan to remove the support of gjk_guess in
  ///        CollisionRequest. If you want to use this, access it by the
  ///        GJKSolver.
  bool is_gjk_guess_valid{false};
  Vector3<S> gjk_guess;

  friend std::ostream& operator<<(std::ostream& out, const GJKSolver& solver) {
    out << "GJKSolver"
        << "\n    gjk tolerance:       " << solver.gjk_tolerance
        << "\n    gjk max iterations:  " << solver.gjk_max_iterations
        << "\n    epa tolerance:       " << solver.epa_tolerance
        << "\n    epa max face num:    " << solver.epa_max_face_num
        << "\n    epa max vertex num:  " << solver.epa_max_vertex_num
        << "\n    epa max iterations:  " << solver.epa_max_iterations
        << "\n    enable cached guess: " << solver.is_gjk_guess_valid;
    if (solver.is_gjk_guess_valid) out << solver.gjk_guess.transpose();
    return out;
  }
};

using GJKSolverf = GJKSolver<float>;
using GJKSolverd = GJKSolver<double>;

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/gjk_solver-inl.h"

#endif
