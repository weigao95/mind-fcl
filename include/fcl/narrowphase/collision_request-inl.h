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

#ifndef FCL_COLLISIONREQUEST_INL_H
#define FCL_COLLISIONREQUEST_INL_H

#include "fcl/narrowphase/collision_request.h"
#include "fcl/narrowphase/collision_result.h"

namespace fcl {

//==============================================================================
template <typename S>
CollisionRequest<S>::CollisionRequest() : CollisionRequest(1) {
  // Do nothing
}

//==============================================================================
template <typename S>
CollisionRequest<S>::CollisionRequest(std::size_t n_max_contacts)
    : num_max_contacts_(n_max_contacts),
      penetration_mode_(detail::CollisionPenetrationType::Disabled),
      binary_collision_tolerance_(Real(1e-6)),
      distance_tolerance_(Real(1e-6)) {
  // Do nothing
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::disablePenetration() {
  penetration_mode_.setAsDisabled();
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::useDefaultPenetration() {
  penetration_mode_.setAsDefaultGJK_EPA();
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::useMinimumDistancePenetration() {
  penetration_mode_.setAsDefaultGJK_EPA();
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::useDirectedPenetration(
    Vector3<S> shape2_escape_direction) {
  const S direction_norm = shape2_escape_direction.norm();
  if (shape2_escape_direction.norm() <= S(0.0)) {
    shape2_escape_direction = Vector3<S>::UnitZ();
  } else {
    shape2_escape_direction /= direction_norm;
  }
  penetration_mode_.setAsDirectedPenetration(
      std::move(shape2_escape_direction));
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::useIncrementalMinimumDistancePenetration(
    Vector3<S> shape2_escape_direction) {
  const S direction_norm = shape2_escape_direction.norm();
  if (shape2_escape_direction.norm() <= S(0.0)) {
    shape2_escape_direction = Vector3<S>::UnitZ();
  } else {
    shape2_escape_direction /= direction_norm;
  }
  penetration_mode_.setAsIncrementalMinimumPenetration(
      std::move(shape2_escape_direction));
}

//==============================================================================
template <typename S>
bool CollisionRequest<S>::isPenetrationEnabled() const {
  return penetration_mode_.penetration_type() !=
         detail::CollisionPenetrationType::Disabled;
}

//==============================================================================
template <typename S>
const detail::CollisionPenetrationMode<S>&
CollisionRequest<S>::penetrationMode() const {
  return penetration_mode_;
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::setMaxContactCount(std::size_t n_max_contact) {
  num_max_contacts_ = n_max_contact;
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::setBinaryCollisionTolerance(Real tolerance) {
  binary_collision_tolerance_ = tolerance;
}

//==============================================================================
template <typename S>
void CollisionRequest<S>::setPenetrationDistanceTolerance(Real tolerance) {
  distance_tolerance_ = tolerance;
}

//==============================================================================
template <typename S>
bool CollisionRequest<S>::terminationConditionSatisfied(
    const CollisionResult<S>& result) const {
  return result.terminationConditionSatisfied(num_max_contacts_);
}

}  // namespace fcl

#endif
