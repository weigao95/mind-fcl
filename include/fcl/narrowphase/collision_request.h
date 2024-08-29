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

#ifndef FCL_COLLISIONREQUEST_H
#define FCL_COLLISIONREQUEST_H

#include "fcl/common/types.h"
#include "fcl/narrowphase/detail/collision_penetration_mode.h"

namespace fcl {

template <typename S>
struct CollisionResult;

/// @brief Parameters for performing collision request.
template <typename S>
struct CollisionRequest {
 private:
  /// @brief The maximum number of contacts that can be returned.
  std::size_t num_max_contacts_;

  /// @brief If true, contact information (e.g., normal, penetration depth, and
  /// contact position) will be returned.
  detail::CollisionPenetrationMode<S> penetration_mode_;

  /// @brief Numerical tolerance to use in the GJK algorithm.
  /// NOTE: The default value is currently set as 1e-6 to provide backwards
  /// compatibility; historically it has been 1e-6. Future code should provide
  /// a value that is consistent with the precision of `S`.
  using Real = typename Eigen::NumTraits<S>::Real;
  Real binary_collision_tolerance_{1e-6};
  Real distance_tolerance_{1e-6};

 public:
  /// @brief Default constructor
  explicit CollisionRequest();
  explicit CollisionRequest(std::size_t n_max_contacts);

  /// @brief Contact type related option
  void disablePenetration();
  void useDefaultPenetration();
  void useMinimumDistancePenetration();
  void useDirectedPenetration(Vector3<S> shape2_escape_direction_in_world);
  void useIncrementalMinimumDistancePenetration(
      Vector3<S> shape2_escape_direction_in_world_hint);
  bool isPenetrationEnabled() const;
  const detail::CollisionPenetrationMode<S>& penetrationMode() const;

  /// @brief Contact count related option
  std::size_t maxNumContacts() const { return num_max_contacts_; }
  void setMaxContactCount(std::size_t n_max_contact);

  /// Tolerance
  void setBinaryCollisionTolerance(Real tolerance);
  void setPenetrationDistanceTolerance(Real tolerance);
  Real binaryCollisionTolerance() const { return binary_collision_tolerance_; };
  Real distanceTolerance() const { return distance_tolerance_; }

  /// Whether the termination condition is meet
  bool terminationConditionSatisfied(const CollisionResult<S>& result) const;
};

using CollisionRequestf = CollisionRequest<float>;
using CollisionRequestd = CollisionRequest<double>;

}  // namespace fcl

#include "fcl/narrowphase/collision_request-inl.h"

#endif
