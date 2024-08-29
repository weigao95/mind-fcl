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

#ifndef FCL_NARROWPHASE_GJKSOLVERINDEP_INL_H
#define FCL_NARROWPHASE_GJKSOLVERINDEP_INL_H

#include <algorithm>

#include "fcl/common/unused.h"
#include "fcl/cvx_collide/epa.h"
#include "fcl/cvx_collide/gjk.h"
#include "fcl/cvx_collide/mpr.h"
#include "fcl/geometry/shape/triangle_p.h"
#include "fcl/narrowphase/detail/gjk_solver.h"
#include "fcl/narrowphase/detail/gjk_solver_cvx.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_triangle.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/halfspace.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/plane.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_capsule.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_cylinder.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_sphere.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_triangle.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/tetrahedron_intersect.h"

namespace fcl {
namespace detail {

//==============================================================================
// Global compute option
constexpr bool use_mpr_if_no_contact = true;

//==============================================================================
template <typename S, typename Shape1, typename Shape2>
struct ShapeIntersectIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape1& s1,
                  const Transform3<S>& tf1, const Shape2& s2,
                  const Transform3<S>& tf2,
                  CollisionPenetrationContactData<S>* contacts) {
    // Construct the shape
    Vector3<S> guess(1, 0, 0);
    if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;
    MinkowskiDiff<S> shape;
    shape.shapes[0] = constructGJKGeometry(&s1);
    shape.shapes[1] = constructGJKGeometry(&s2);
    shape.support_function = computeSupport<S>;
    shape.interior_function = computeInterior<S>;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    // Use mpr if no contact required
    if (use_mpr_if_no_contact && (contacts == nullptr)) {
      MPR<S> mpr(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
      auto mpr_result = mpr.Intersect(shape);

      // Check the result
      if (mpr_result == MPR<S>::IntersectStatus::Intersect) {
        return true;
      } else if (mpr_result == MPR<S>::IntersectStatus::Separated) {
        return false;
      }  // else failed (very low probability), but move to gjk
    }

    // Check can we use GJK only
    GJK2<S> gjk2(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    GJKSimplex<S> gjk_simplex;
    auto gjk_result = gjk2.Evaluate(shape, gjk_simplex, -guess);
    if (contacts == nullptr) {
      return gjk_result == GJK_Status::Intersect;
    }

    // We must run EPA for contact
    assert(contacts != nullptr);
    if (gjk_result == GJK_Status::Intersect) {
      EPA2<S> epa2(gjk_solver.epa_max_face_num, gjk_solver.epa_max_iterations,
                   gjk_solver.epa_tolerance);
      S epa_depth{0};
      Vector3<S> p_on_1, p_on_2;
      auto epa_status_out =
          epa2.Evaluate(gjk_simplex, shape, &epa_depth, &p_on_1, &p_on_2);
      if (epa_status_out != EPA_Status::Failed) {
        Vector3<S> normal_in_1 = p_on_1 - p_on_2;
        if (normal_in_1.squaredNorm() <= 0)
          normal_in_1 = Vector3<S>::UnitZ();
        else
          normal_in_1.normalize();
        Vector3<S> point_in_1 = S(0.5) * (p_on_1 + p_on_2);
        Vector3<S> point = tf1 * point_in_1;
        Vector3<S> normal = tf1.linear().matrix() * normal_in_1;
        // Note: use position distance
        const S depth = epa_depth;
        contacts->emplace_back(normal, point, depth);
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
};

//==============================================================================
template <typename S>
template <typename Shape1, typename Shape2>
bool GJKSolver<S>::shapeIntersect(
    const Shape1& s1, const Transform3<S>& tf1, const Shape2& s2,
    const Transform3<S>& tf2,
    CollisionPenetrationContactData<S>* contacts) const {
  return ShapeIntersectIndepImpl<S, Shape1, Shape2>::run(*this, s1, tf1, s2,
                                                         tf2, contacts);
}

// clang-format off
// Shape intersect algorithms not using built-in GJK algorithm
//
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// |            | box | sphere | ellipsoid | capsule | cone | cylinder | plane | half-space | triangle |  convex  |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | box        |  O  |   O    |           |         |      |          |   O   |      O     |          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | sphere     |/////|   O    |           |    O    |      |    O     |   O   |      O     |     O    |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | ellipsoid  |/////|////////|           |         |      |          |   O   |      O     |          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | capsule    |/////|////////|///////////|         |      |          |   O   |      O     |          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | cone       |/////|////////|///////////|/////////|      |          |   O   |      O     |          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | cylinder   |/////|////////|///////////|/////////|//////|          |   O   |      O     |          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | plane      |/////|////////|///////////|/////////|//////|//////////|   O   |      O     |     O    |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | half-space |/////|////////|///////////|/////////|//////|//////////|///////|      O     |     O    |    O     |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | triangle   |/////|////////|///////////|/////////|//////|//////////|///////|////////////|          |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// | convex     |/////|////////|///////////|/////////|//////|//////////|///////|////////////|//////////|          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+----------+
// clang-format on

#define FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT_REG(SHAPE1, SHAPE2, ALG)         \
  template <typename S>                                                      \
  struct ShapeIntersectIndepImpl<S, SHAPE1<S>, SHAPE2<S>> {                  \
    static bool run(const GJKSolver<S>& /*gjk_solver*/, const SHAPE1<S>& s1, \
                    const Transform3<S>& tf1, const SHAPE2<S>& s2,           \
                    const Transform3<S>& tf2,                                \
                    CollisionPenetrationContactData<S>* contacts) {          \
      return ALG(s1, tf1, s2, tf2, contacts);                                \
    }                                                                        \
  };

#define FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT_INV(SHAPE1, SHAPE2, ALG)         \
  template <typename S>                                                      \
  struct ShapeIntersectIndepImpl<S, SHAPE2<S>, SHAPE1<S>> {                  \
    static bool run(const GJKSolver<S>& /*gjk_solver*/, const SHAPE2<S>& s1, \
                    const Transform3<S>& tf1, const SHAPE1<S>& s2,           \
                    const Transform3<S>& tf2,                                \
                    CollisionPenetrationContactData<S>* contacts) {          \
      const bool res = ALG(s2, tf2, s1, tf1, contacts);                      \
      if (contacts) flipNormal(*contacts);                                   \
      return res;                                                            \
    }                                                                        \
  };

#define FCL_GJK_INDEP_SHAPE_INTERSECT(SHAPE, ALG) \
  FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT_REG(SHAPE, SHAPE, ALG)

#define FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(SHAPE1, SHAPE2, ALG) \
  FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT_REG(SHAPE1, SHAPE2, ALG)   \
  FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT_INV(SHAPE1, SHAPE2, ALG)

FCL_GJK_INDEP_SHAPE_INTERSECT(Sphere, detail::sphereSphereIntersect)
FCL_GJK_INDEP_SHAPE_INTERSECT(Box, detail::boxBoxIntersect)

FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Sphere, Capsule,
                                    detail::sphereCapsuleIntersect)

FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Sphere, Box, detail::sphereBoxIntersect)

FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Sphere, Cylinder,
                                    detail::sphereCylinderIntersect)

FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Sphere, Halfspace,
                                    detail::sphereHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Ellipsoid, Halfspace,
                                    detail::ellipsoidHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Box, Halfspace,
                                    detail::boxHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Capsule, Halfspace,
                                    detail::capsuleHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Cylinder, Halfspace,
                                    detail::cylinderHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Cone, Halfspace,
                                    detail::coneHalfspaceIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Convex, Halfspace,
                                    detail::convexHalfspaceIntersect)

FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Sphere, Plane, detail::spherePlaneIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Ellipsoid, Plane,
                                    detail::ellipsoidPlaneIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Box, Plane, detail::boxPlaneIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Capsule, Plane,
                                    detail::capsulePlaneIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Cylinder, Plane,
                                    detail::cylinderPlaneIntersect)
FCL_GJK_INDEP_SHAPE_SHAPE_INTERSECT(Cone, Plane, detail::conePlaneIntersect)

template <typename S>
struct ShapeIntersectIndepImpl<S, Halfspace<S>, Halfspace<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Halfspace<S>& s1,
                  const Transform3<S>& tf1, const Halfspace<S>& s2,
                  const Transform3<S>& tf2,
                  CollisionPenetrationContactData<S>* contacts) {
    FCL_UNUSED(contacts);

    Halfspace<S> s;
    Vector3<S> p, d;
    S depth;
    int ret;
    return detail::halfspaceIntersect(s1, tf1, s2, tf2, p, d, s, depth, ret);
  }
};

template <typename S>
struct ShapeIntersectIndepImpl<S, Plane<S>, Plane<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Plane<S>& s1,
                  const Transform3<S>& tf1, const Plane<S>& s2,
                  const Transform3<S>& tf2,
                  CollisionPenetrationContactData<S>* contacts) {
    return detail::planeIntersect(s1, tf1, s2, tf2, contacts);
  }
};

template <typename S>
struct ShapeIntersectIndepImpl<S, Plane<S>, Halfspace<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Plane<S>& s1,
                  const Transform3<S>& tf1, const Halfspace<S>& s2,
                  const Transform3<S>& tf2,
                  CollisionPenetrationContactData<S>* contacts) {
    FCL_UNUSED(contacts);

    Plane<S> pl;
    Vector3<S> p, d;
    S depth;
    int ret;
    return detail::planeHalfspaceIntersect(s1, tf1, s2, tf2, pl, p, d, depth,
                                           ret);
  }
};

template <typename S>
struct ShapeIntersectIndepImpl<S, Halfspace<S>, Plane<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Halfspace<S>& s1,
                  const Transform3<S>& tf1, const Plane<S>& s2,
                  const Transform3<S>& tf2,
                  CollisionPenetrationContactData<S>* contacts) {
    FCL_UNUSED(contacts);

    Plane<S> pl;
    Vector3<S> p, d;
    S depth;
    int ret;
    return detail::halfspacePlaneIntersect(s1, tf1, s2, tf2, pl, p, d, depth,
                                           ret);
  }
};

//==============================================================================
template <typename S, typename Shape>
bool shapeTriangleIntersectImplRun(const GJKSolver<S>& gjk_solver,
                                   const Shape& s, const Transform3<S>& tf,
                                   const Vector3<S>& P1, const Vector3<S>& P2,
                                   const Vector3<S>& P3,
                                   Vector3<S>* contact_points,
                                   S* penetration_depth, Vector3<S>* normal) {
  // Construct the triangle
  TriangleP<S> tri(P1, P2, P3);
  Vector3<S> guess(1, 0, 0);
  if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

  // Make the space
  MinkowskiDiff<S> shape;
  shape.shapes[0] = constructGJKGeometry(&s);
  shape.shapes[1] = constructGJKGeometry(&tri);
  shape.support_function = computeSupport<S>;
  shape.interior_function = computeInterior<S>;
  shape.toshape1 = tf.linear();
  shape.toshape0 = tf.inverse(Eigen::Isometry);

  // The compute option allows MPR
  if (use_mpr_if_no_contact && contact_points == nullptr &&
      penetration_depth == nullptr && normal == nullptr) {
    MPR<S> mpr(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    auto mpr_result = mpr.Intersect(shape);
    return mpr_result == MPR<S>::IntersectStatus::Intersect;
  }

  // Check can we use GJK only
  GJK2<S> gjk2(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
  GJKSimplex<S> gjk_simplex;
  auto gjk_result = gjk2.Evaluate(shape, gjk_simplex, -guess);
  if (contact_points == nullptr && penetration_depth == nullptr &&
      normal == nullptr) {
    return gjk_result == GJK_Status::Intersect;
  }

  // We must run EPA
  if (gjk_result == GJK_Status::Intersect) {
    EPA2<S> epa2(gjk_solver.epa_max_face_num, gjk_solver.epa_max_iterations,
                 gjk_solver.epa_tolerance);
    S epa_depth{0};
    Vector3<S> p_on_0, p_on_1;
    auto epa_status_out =
        epa2.Evaluate(gjk_simplex, shape, &epa_depth, &p_on_0, &p_on_1);
    if (epa_status_out != EPA_Status::Failed) {
      // Assign the contact point, the middle of two points
      if (contact_points != nullptr) {
        Vector3<S> point_in_0 = S(0.5) * (p_on_0 + p_on_1);
        *contact_points = tf * point_in_0;
      }

      // The normal in world
      if (normal != nullptr) {
        Vector3<S> normal_in_0 = p_on_0 - p_on_1;
        if (normal_in_0.squaredNorm() <= 0)
          normal_in_0 = Vector3<S>::UnitZ();
        else
          normal_in_0.normalize();
        *normal = tf.linear().matrix() * normal_in_0;
      }

      // The penetration depth
      if (penetration_depth != nullptr) {
        *penetration_depth = epa_depth;
      }

      // Done
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

//==============================================================================
template <typename S, typename Shape>
struct ShapeTriangleIntersectIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape& s,
                  const Transform3<S>& tf, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  Vector3<S>* contact_points, S* penetration_depth,
                  Vector3<S>* normal) {
    return shapeTriangleIntersectImplRun<S, Shape>(gjk_solver, s, tf, P1, P2,
                                                   P3, contact_points,
                                                   penetration_depth, normal);
  }
};

template <typename S>
template <typename Shape>
bool GJKSolver<S>::shapeTriangleIntersect(
    const Shape& s, const Transform3<S>& tf, const Vector3<S>& P1,
    const Vector3<S>& P2, const Vector3<S>& P3,
    CollisionPenetrationContactData<S>* penetration_contacts) const {
  if (penetration_contacts == nullptr) {
    return ShapeTriangleIntersectIndepImpl<S, Shape>::run(
        *this, s, tf, P1, P2, P3, nullptr, nullptr, nullptr);
  } else {
    ContactPoint<S> contact;
    const auto intersect = ShapeTriangleIntersectIndepImpl<S, Shape>::run(
        *this, s, tf, P1, P2, P3, &contact.pos, &contact.penetration_depth,
        &contact.normal);
    if (intersect) {
      penetration_contacts->emplace_back(std::move(contact));
    }
    return intersect;
  }
}

//==============================================================================
template <typename S>
struct ShapeTriangleIntersectIndepImpl<S, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s,
                  const Transform3<S>& tf, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  Vector3<S>* contact_points, S* penetration_depth,
                  Vector3<S>* normal) {
    return detail::sphereTriangleIntersect(s, tf, P1, P2, P3, contact_points,
                                           penetration_depth, normal);
  }
};

//==============================================================================
template <typename S>
struct ShapeTriangleIntersectIndepImpl<S, Box<S>> {
  static bool run(const GJKSolver<S>& gjk_solver, const Box<S>& s,
                  const Transform3<S>& tf, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  Vector3<S>* contact_points, S* penetration_depth,
                  Vector3<S>* normal) {
    if (contact_points == nullptr && penetration_depth == nullptr &&
        normal == nullptr) {
      return detail::boxTriangleIntersect(s, tf, P1, P2, P3);
    }
    return shapeTriangleIntersectImplRun<S>(gjk_solver, s, tf, P1, P2, P3,
                                            contact_points, penetration_depth,
                                            normal);
  }
};

//==============================================================================
template <typename S, typename Shape>
bool shapeTransformedTriangleIntersectIndepImplRun(
    const GJKSolver<S>& gjk_solver, const Shape& s, const Transform3<S>& tf1,
    const Vector3<S>& P1, const Vector3<S>& P2, const Vector3<S>& P3,
    const Transform3<S>& tf2, Vector3<S>* contact_points, S* penetration_depth,
    Vector3<S>* normal) {
  // Construct the triangle
  TriangleP<S> tri(P1, P2, P3);
  Vector3<S> guess(1, 0, 0);
  if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

  // The shape and transform
  MinkowskiDiff<S> shape;
  shape.shapes[0] = constructGJKGeometry(&s);
  shape.shapes[1] = constructGJKGeometry(&tri);
  shape.support_function = computeSupport<S>;
  shape.interior_function = computeInterior<S>;
  shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
  shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

  // The compute option allows MPR only
  if (use_mpr_if_no_contact && contact_points == nullptr &&
      penetration_depth == nullptr && normal == nullptr) {
    MPR<S> mpr(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    auto mpr_result = mpr.Intersect(shape);
    return mpr_result == MPR<S>::IntersectStatus::Intersect;
  }

  // Check can we use GJK only
  GJK2<S> gjk2(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
  GJKSimplex<S> gjk_simplex;
  auto gjk_result = gjk2.Evaluate(shape, gjk_simplex, -guess);
  if (contact_points == nullptr && penetration_depth == nullptr &&
      normal == nullptr) {
    return gjk_result == GJK_Status::Intersect;
  }

  // We must run EPA
  if (gjk_result == GJK_Status::Intersect) {
    EPA2<S> epa2(gjk_solver.epa_max_face_num, gjk_solver.epa_max_iterations,
                 gjk_solver.epa_tolerance);
    S epa_depth{0};
    Vector3<S> p_on_1, p_on_2;
    auto epa_status_out =
        epa2.Evaluate(gjk_simplex, shape, &epa_depth, &p_on_1, &p_on_2);
    if (epa_status_out != EPA_Status::Failed) {
      // Assign the contact point, the middle of two points
      if (contact_points != nullptr) {
        Vector3<S> point_in_1 = S(0.5) * (p_on_1 + p_on_2);
        *contact_points = tf1 * point_in_1;
      }

      // The normal in world
      if (normal != nullptr) {
        Vector3<S> normal_in_1 = p_on_1 - p_on_2;
        if (normal_in_1.squaredNorm() <= 0)
          normal_in_1 = Vector3<S>::UnitZ();
        else
          normal_in_1.normalize();
        *normal = tf1.linear().matrix() * normal_in_1;
      }

      // The penetration depth
      if (penetration_depth != nullptr) {
        *penetration_depth = epa_depth;
      }

      // Done
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

//==============================================================================
template <typename S, typename Shape>
struct ShapeTransformedTriangleIntersectIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, Vector3<S>* contact_points,
                  S* penetration_depth, Vector3<S>* normal) {
    return shapeTransformedTriangleIntersectIndepImplRun<S, Shape>(
        gjk_solver, s, tf1, P1, P2, P3, tf2, contact_points, penetration_depth,
        normal);
  }
};

//==============================================================================
template <typename S>
template <typename Shape>
bool GJKSolver<S>::shapeTriangleIntersect(
    const Shape& s, const Transform3<S>& tf1, const Vector3<S>& P1,
    const Vector3<S>& P2, const Vector3<S>& P3, const Transform3<S>& tf2,
    CollisionPenetrationContactData<S>* penetration_contacts) const {
  if (penetration_contacts == nullptr) {
    return ShapeTransformedTriangleIntersectIndepImpl<S, Shape>::run(
        *this, s, tf1, P1, P2, P3, tf2, nullptr, nullptr, nullptr);
  } else {
    ContactPoint<S> contact;
    const bool intersect =
        ShapeTransformedTriangleIntersectIndepImpl<S, Shape>::run(
            *this, s, tf1, P1, P2, P3, tf2, &contact.pos,
            &contact.penetration_depth, &contact.normal);
    if (intersect) {
      penetration_contacts->emplace_back(std::move(contact));
    }
    return intersect;
  }
}

//==============================================================================
template <typename S>
struct ShapeTransformedTriangleIntersectIndepImpl<S, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, Vector3<S>* contact_points,
                  S* penetration_depth, Vector3<S>* normal) {
    return detail::sphereTriangleIntersect(s, tf1, tf2 * P1, tf2 * P2, tf2 * P3,
                                           contact_points, penetration_depth,
                                           normal);
  }
};

//==============================================================================
template <typename S>
struct ShapeTransformedTriangleIntersectIndepImpl<S, Box<S>> {
  static bool run(const GJKSolver<S>& gjk_solver, const Box<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, Vector3<S>* contact_points,
                  S* penetration_depth, Vector3<S>* normal) {
    if (contact_points == nullptr && penetration_depth == nullptr &&
        normal == nullptr) {
      return detail::boxTriangleIntersect(s, tf1, P1, P2, P3, tf2);
    }
    return shapeTransformedTriangleIntersectIndepImplRun<S>(
        gjk_solver, s, tf1, P1, P2, P3, tf2, contact_points, penetration_depth,
        normal);
  }
};

//==============================================================================
template <typename S>
struct ShapeTransformedTriangleIntersectIndepImpl<S, Halfspace<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Halfspace<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, Vector3<S>* contact_points,
                  S* penetration_depth, Vector3<S>* normal) {
    return detail::halfspaceTriangleIntersect(
        s, tf1, P1, P2, P3, tf2, contact_points, penetration_depth, normal);
  }
};

//==============================================================================
template <typename S>
struct ShapeTransformedTriangleIntersectIndepImpl<S, Plane<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Plane<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, Vector3<S>* contact_points,
                  S* penetration_depth, Vector3<S>* normal) {
    return detail::planeTriangleIntersect(
        s, tf1, P1, P2, P3, tf2, contact_points, penetration_depth, normal);
  }
};

//==============================================================================
template <typename S, typename Shape>
bool shapeTransformedTetrahedronIntersectIndepImplRun(
    const GJKSolver<S>& gjk_solver, const Shape& s, const Transform3<S>& tf1,
    const Vector3<S>& P1, const Vector3<S>& P2, const Vector3<S>& P3,
    const Vector3<S>& P4, const Transform3<S>& tf2, Vector3<S>* contact_points,
    S* penetration_depth, Vector3<S>* normal) {
  // Construct the triangle
  Tetrahedron<S> tet(P1, P2, P3, P4);
  Vector3<S> guess(1, 0, 0);
  if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

  // The shape and transform
  MinkowskiDiff<S> shape;
  shape.shapes[0] = constructGJKGeometry(&s);
  shape.shapes[1] = constructGJKGeometry(&tet);
  shape.support_function = computeSupport<S>;
  shape.interior_function = computeInterior<S>;
  shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
  shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

  // The compute option allows MPR only
  if (use_mpr_if_no_contact && contact_points == nullptr &&
      penetration_depth == nullptr && normal == nullptr) {
    MPR<S> mpr(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    auto mpr_result = mpr.Intersect(shape);
    return mpr_result == MPR<S>::IntersectStatus::Intersect;
  }

  // Check can we use GJK only
  GJK2<S> gjk2(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
  GJKSimplex<S> gjk_simplex;
  auto gjk_result = gjk2.Evaluate(shape, gjk_simplex, -guess);
  if (contact_points == nullptr && penetration_depth == nullptr &&
      normal == nullptr) {
    return gjk_result == GJK_Status::Intersect;
  }

  // We must run EPA
  if (gjk_result == GJK_Status::Intersect) {
    EPA2<S> epa2(gjk_solver.epa_max_face_num, gjk_solver.epa_max_iterations,
                 gjk_solver.epa_tolerance);
    S epa_depth{0};
    Vector3<S> p_on_1, p_on_2;
    auto epa_status_out =
        epa2.Evaluate(gjk_simplex, shape, &epa_depth, &p_on_1, &p_on_2);
    if (epa_status_out != EPA_Status::Failed) {
      // Assign the contact point, the middle of two points
      if (contact_points != nullptr) {
        Vector3<S> point_in_1 = S(0.5) * (p_on_1 + p_on_2);
        *contact_points = tf1 * point_in_1;
      }

      // The normal in world
      if (normal != nullptr) {
        Vector3<S> normal_in_1 = p_on_1 - p_on_2;
        if (normal_in_1.squaredNorm() <= 0)
          normal_in_1 = Vector3<S>::UnitZ();
        else
          normal_in_1.normalize();
        *normal = tf1.linear().matrix() * normal_in_1;
      }

      // The penetration depth
      if (penetration_depth != nullptr) {
        *penetration_depth = epa_depth;
      }

      // Done
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

//==============================================================================
template <typename S, typename Shape>
struct ShapeTransformedTetrahedronIntersectIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Vector3<S>& P4, const Transform3<S>& tf2,
                  Vector3<S>* contact_points, S* penetration_depth,
                  Vector3<S>* normal) {
    return shapeTransformedTetrahedronIntersectIndepImplRun<S, Shape>(
        gjk_solver, s, tf1, P1, P2, P3, P4, tf2, contact_points,
        penetration_depth, normal);
  }
};

//==============================================================================
template <typename S>
struct ShapeTransformedTetrahedronIntersectIndepImpl<S, Box<S>> {
  static bool run(const GJKSolver<S>& gjk_solver, const Box<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Vector3<S>& P4, const Transform3<S>& tf2,
                  Vector3<S>* contact_points, S* penetration_depth,
                  Vector3<S>* normal) {
    if (contact_points == nullptr && penetration_depth == nullptr &&
        normal == nullptr) {
      Tetrahedron<S> tetrahedron(P1, P2, P3, P4);
      return detail::boxTerahedronIntersect(s, tf1, tetrahedron, tf2);
    }
    return shapeTransformedTetrahedronIntersectIndepImplRun<S>(
        gjk_solver, s, tf1, P1, P2, P3, P4, tf2, contact_points,
        penetration_depth, normal);
  }
};

//==============================================================================
template <typename S>
template <typename Shape>
bool GJKSolver<S>::shapeTetrahedronIntersect(
    const Shape& s, const Transform3<S>& tf1, const Vector3<S>& P1,
    const Vector3<S>& P2, const Vector3<S>& P3, const Vector3<S>& P4,
    const Transform3<S>& tf2,
    CollisionPenetrationContactData<S>* penetration_contacts) const {
  if (penetration_contacts == nullptr) {
    return ShapeTransformedTetrahedronIntersectIndepImpl<S, Shape>::run(
        *this, s, tf1, P1, P2, P3, P4, tf2, nullptr, nullptr, nullptr);
  } else {
    ContactPoint<S> contact;
    const bool intersect =
        ShapeTransformedTetrahedronIntersectIndepImpl<S, Shape>::run(
            *this, s, tf1, P1, P2, P3, P4, tf2, &contact.pos,
            &contact.penetration_depth, &contact.normal);
    if (intersect) {
      penetration_contacts->emplace_back(std::move(contact));
    }
    return intersect;
  }
}

//==============================================================================
template <typename S, typename Shape1, typename Shape2>
struct ShapeDistanceIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape1& s1,
                  const Transform3<S>& tf1, const Shape2& s2,
                  const Transform3<S>& tf2, S* distance, Vector3<S>* p1,
                  Vector3<S>* p2) {
    Vector3<S> guess(1, 0, 0);
    if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

    // Construct the shape
    MinkowskiDiff<S> shape;
    shape.shapes[0] = constructGJKGeometry(&s1);
    shape.shapes[1] = constructGJKGeometry(&s2);
    shape.support_function = computeSupport<S>;
    shape.interior_function = computeInterior<S>;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    // Run gjk
    GJK2<S> gjk(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto gjk_result =
        gjk.Evaluate(shape, simplex, -guess, &min_distance_output);
    if (gjk_result == GJK_Status::Separated) {
      // Answer is solved in Shape1's local frame; answers are given in the
      // world frame.
      if (p1) p1->noalias() = tf1 * min_distance_output.p0_if_separated;
      if (p2) p2->noalias() = tf1 * min_distance_output.p1_if_separated;
      if (distance) *distance = min_distance_output.separation_distance();
      return true;
    } else {
      if (distance) *distance = -1;
      return false;
    }
  }
};

template <typename S>
template <typename Shape1, typename Shape2>
bool GJKSolver<S>::shapeDistance(const Shape1& s1, const Transform3<S>& tf1,
                                 const Shape2& s2, const Transform3<S>& tf2,
                                 S* dist, Vector3<S>* p1,
                                 Vector3<S>* p2) const {
  return ShapeDistanceIndepImpl<S, Shape1, Shape2>::run(*this, s1, tf1, s2, tf2,
                                                        dist, p1, p2);
}

//==============================================================================
template <typename S, typename Shape1, typename Shape2>
struct ShapeSignedDistanceIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape1& s1,
                  const Transform3<S>& tf1, const Shape2& s2,
                  const Transform3<S>& tf2, S* distance, Vector3<S>* p1,
                  Vector3<S>* p2) {
    Vector3<S> guess(1, 0, 0);
    if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

    // Construct the shape
    MinkowskiDiff<S> shape;
    shape.shapes[0] = constructGJKGeometry(&s1);
    shape.shapes[1] = constructGJKGeometry(&s2);
    shape.support_function = computeSupport<S>;
    shape.interior_function = computeInterior<S>;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    // Run gjk
    GJK2<S> gjk(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto gjk_result =
        gjk.Evaluate(shape, simplex, -guess, &min_distance_output);
    if (gjk_result == GJK_Status::Separated &&
        (min_distance_output.is_separation_point_valid)) {
      // Answer is solved in Shape1's local frame; answers are given in the
      // world frame.
      if (p1) p1->noalias() = tf1 * min_distance_output.p0_if_separated;
      if (p2) p2->noalias() = tf1 * min_distance_output.p1_if_separated;
      if (distance) *distance = min_distance_output.separation_distance();
      return true;
    } else if (gjk_result == GJK_Status::Intersect) {
      // Invoke epa
      EPA2<S> epa2(gjk_solver.epa_max_face_num, gjk_solver.epa_max_iterations,
                   gjk_solver.epa_tolerance);
      Vector3<S> p_on_1, p_on_2;
      S depth_if_penetration;
      auto epa_status = epa2.Evaluate(simplex, shape, &depth_if_penetration,
                                      &p_on_1, &p_on_2);
      if (epa_status != EPA_Status::Failed) {
        if (p1 != nullptr) *p1 = tf1 * p_on_1;
        if (p2 != nullptr) *p2 = tf1 * p_on_2;
        if (distance != nullptr) *distance = -depth_if_penetration;
        return true;
      } else {
        return false;
      }
    } else {
      if (distance) *distance = -1;
      return false;
    }
  }
};

template <typename S>
template <typename Shape1, typename Shape2>
bool GJKSolver<S>::shapeSignedDistance(const Shape1& s1,
                                       const Transform3<S>& tf1,
                                       const Shape2& s2,
                                       const Transform3<S>& tf2, S* dist,
                                       Vector3<S>* p1, Vector3<S>* p2) const {
  return ShapeSignedDistanceIndepImpl<S, Shape1, Shape2>::run(
      *this, s1, tf1, s2, tf2, dist, p1, p2);
}

// Shape distance algorithms not using built-in GJK algorithm
// clang-format off
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// |            | box | sphere | ellipsoid | capsule | cone | cylinder | plane | half-space | triangle |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | box        |     |   O    |           |         |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | sphere     |/////|   O    |           |    O    |      |    O     |       |            |     O    |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | ellipsoid  |/////|////////|           |         |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | capsule    |/////|////////|///////////|    O    |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | cone       |/////|////////|///////////|/////////|      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | cylinder   |/////|////////|///////////|/////////|//////|          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | plane      |/////|////////|///////////|/////////|//////|//////////|       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | half-space |/////|////////|///////////|/////////|//////|//////////|///////|            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | triangle   |/////|////////|///////////|/////////|//////|//////////|///////|////////////|          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// clang-format on
//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Sphere<S>, Box<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s1,
                  const Transform3<S>& tf1, const Box<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereBoxDistance(s1, tf1, s2, tf2, dist, p1, p2);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Box<S>, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Box<S>& s1,
                  const Transform3<S>& tf1, const Sphere<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereBoxDistance(s2, tf2, s1, tf1, dist, p2, p1);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Sphere<S>, Capsule<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s1,
                  const Transform3<S>& tf1, const Capsule<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereCapsuleDistance(s1, tf1, s2, tf2, dist, p1, p2);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Capsule<S>, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Capsule<S>& s1,
                  const Transform3<S>& tf1, const Sphere<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereCapsuleDistance(s2, tf2, s1, tf1, dist, p2, p1);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Sphere<S>, Cylinder<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s1,
                  const Transform3<S>& tf1, const Cylinder<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereCylinderDistance(s1, tf1, s2, tf2, dist, p1, p2);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Cylinder<S>, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Cylinder<S>& s1,
                  const Transform3<S>& tf1, const Sphere<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereCylinderDistance(s2, tf2, s1, tf1, dist, p2, p1);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Sphere<S>, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s1,
                  const Transform3<S>& tf1, const Sphere<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereSphereDistance(s1, tf1, s2, tf2, dist, p1, p2);
  }
};

//==============================================================================
template <typename S>
struct ShapeDistanceIndepImpl<S, Capsule<S>, Capsule<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Capsule<S>& s1,
                  const Transform3<S>& tf1, const Capsule<S>& s2,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::capsuleCapsuleDistance(s1, tf1, s2, tf2, dist, p1, p2);
  }
};

//==============================================================================
template <typename S, typename Shape>
struct ShapeTriangleDistanceIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape& s,
                  const Transform3<S>& tf, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3, S* distance,
                  Vector3<S>* p1, Vector3<S>* p2) {
    TriangleP<S> tri(P1, P2, P3);
    Vector3<S> guess(1, 0, 0);
    if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

    // Construct the shape
    MinkowskiDiff<S> shape;
    shape.shapes[0] = constructGJKGeometry(&s);
    shape.shapes[1] = constructGJKGeometry(&tri);
    shape.support_function = computeSupport<S>;
    shape.interior_function = computeInterior<S>;
    shape.toshape1 = tf.linear();
    shape.toshape0 = tf.inverse(Eigen::Isometry);

    // Run gjk
    GJK2<S> gjk(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto gjk_result =
        gjk.Evaluate(shape, simplex, -guess, &min_distance_output);
    if (gjk_result == GJK_Status::Separated &&
        (min_distance_output.is_separation_point_valid)) {
      // Answer is solved in Shape1's local frame; answers are given in the
      // world frame.
      if (p1) p1->noalias() = tf * min_distance_output.p0_if_separated;
      if (p2) p2->noalias() = tf * min_distance_output.p1_if_separated;
      if (distance) *distance = min_distance_output.separation_distance();
      return true;
    } else {
      if (distance) *distance = -1;
      return false;
    }
  }
};

//==============================================================================
template <typename S>
template <typename Shape>
bool GJKSolver<S>::shapeTriangleDistance(const Shape& s,
                                         const Transform3<S>& tf,
                                         const Vector3<S>& P1,
                                         const Vector3<S>& P2,
                                         const Vector3<S>& P3, S* dist,
                                         Vector3<S>* p1, Vector3<S>* p2) const {
  return ShapeTriangleDistanceIndepImpl<S, Shape>::run(*this, s, tf, P1, P2, P3,
                                                       dist, p1, p2);
}

//==============================================================================
template <typename S>
struct ShapeTriangleDistanceIndepImpl<S, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s,
                  const Transform3<S>& tf, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3, S* dist,
                  Vector3<S>* p1, Vector3<S>* p2) {
    return detail::sphereTriangleDistance(s, tf, P1, P2, P3, dist, p1, p2);
  }
};

//==============================================================================
template <typename S, typename Shape>
struct ShapeTransformedTriangleDistanceIndepImpl {
  static bool run(const GJKSolver<S>& gjk_solver, const Shape& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, S* distance, Vector3<S>* p1,
                  Vector3<S>* p2) {
    TriangleP<S> tri(P1, P2, P3);
    Vector3<S> guess(1, 0, 0);
    if (gjk_solver.is_gjk_guess_valid) guess = gjk_solver.gjk_guess;

    // Construct the shape
    MinkowskiDiff<S> shape;
    shape.shapes[0] = constructGJKGeometry(&s);
    shape.shapes[1] = constructGJKGeometry(&tri);
    shape.support_function = computeSupport<S>;
    shape.interior_function = computeInterior<S>;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    // Run gjk
    GJK2<S> gjk(gjk_solver.gjk_max_iterations, gjk_solver.gjk_tolerance);
    GJKSimplex<S> simplex;
    typename GJK2<S>::MinSeparationDistanceOutput min_distance_output;
    auto gjk_result =
        gjk.Evaluate(shape, simplex, -guess, &min_distance_output);
    if (gjk_result == GJK_Status::Separated &&
        (min_distance_output.is_separation_point_valid)) {
      // Answer is solved in Shape1's local frame; answers are given in the
      // world frame.
      if (p1) p1->noalias() = tf1 * min_distance_output.p0_if_separated;
      if (p2) p2->noalias() = tf1 * min_distance_output.p1_if_separated;
      if (distance) *distance = min_distance_output.separation_distance();
      return true;
    } else {
      if (distance) *distance = -1;
      return false;
    }
  }
};

//==============================================================================
template <typename S>
template <typename Shape>
bool GJKSolver<S>::shapeTriangleDistance(
    const Shape& s, const Transform3<S>& tf1, const Vector3<S>& P1,
    const Vector3<S>& P2, const Vector3<S>& P3, const Transform3<S>& tf2,
    S* dist, Vector3<S>* p1, Vector3<S>* p2) const {
  return ShapeTransformedTriangleDistanceIndepImpl<S, Shape>::run(
      *this, s, tf1, P1, P2, P3, tf2, dist, p1, p2);
}

//==============================================================================
template <typename S>
struct ShapeTransformedTriangleDistanceIndepImpl<S, Sphere<S>> {
  static bool run(const GJKSolver<S>& /*gjk_solver*/, const Sphere<S>& s,
                  const Transform3<S>& tf1, const Vector3<S>& P1,
                  const Vector3<S>& P2, const Vector3<S>& P3,
                  const Transform3<S>& tf2, S* dist, Vector3<S>* p1,
                  Vector3<S>* p2) {
    return detail::sphereTriangleDistance(s, tf1, P1, P2, P3, tf2, dist, p1,
                                          p2);
  }
};

//==============================================================================
template <typename S>
GJKSolver<S>::GJKSolver() {
  gjk_max_iterations = 128;
  gjk_tolerance = constants<S>::gjk_default_tolerance();
  epa_max_face_num = 256;
  epa_max_vertex_num = 64;
  epa_max_iterations = 255;
  epa_tolerance = constants<S>::gjk_default_tolerance();
  is_gjk_guess_valid = false;
  gjk_guess = Vector3<S>(1, 0, 0);
}

}  // namespace detail
}  // namespace fcl

#endif
