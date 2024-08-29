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

#ifndef FCL_COLLISIONRESULT_INL_H
#define FCL_COLLISIONRESULT_INL_H

#include "fcl/narrowphase/collision_result.h"

namespace fcl {

//==============================================================================
template <typename S>
CollisionResult<S>::CollisionResult() : user_stop_(false) {
  // Do nothing
}

//==============================================================================
template <typename S>
CollisionResult<S>::CollisionResult(UserContactProcessFunctor<S> user_functor)
    : user_process_functor_(std::move(user_functor)), user_stop_(false) {
  // Do nothing
}

//==============================================================================
template <typename S>
void CollisionResult<S>::addContact(const Contact<S>& c) {
  if (user_process_functor_) {
    bool suggest_stop = false;
    bool keep_this = true;
    user_process_functor_(c, keep_this, suggest_stop);
    user_stop_ |= suggest_stop;
    if (keep_this) {
      contacts_.push_back(c);
    }
  } else {
    contacts_.push_back(c);
  }
}

//==============================================================================
template <typename S>
void CollisionResult<S>::clear() {
  user_stop_ = false;
  contacts_.clear();
}

//==============================================================================
template <typename S>
std::size_t CollisionResult<S>::numContacts() const {
  return contacts_.size();
}

//==============================================================================
template <typename S>
bool CollisionResult<S>::isCollision() const {
  return numContacts() > 0;
}

//==============================================================================
template <typename S>
bool CollisionResult<S>::terminationConditionSatisfied(
    std::size_t n_max_contacts) const {
  return user_stop_ || (isCollision() && numContacts() >= n_max_contacts);
}

//==============================================================================
template <typename S>
const Contact<S>& CollisionResult<S>::getContact(std::size_t i) const {
  if (i < contacts_.size())
    return contacts_[i];
  else
    return contacts_.back();
}

//==============================================================================
template <typename S>
void CollisionResult<S>::getContacts(std::vector<Contact<S>>& contacts_in) {
  contacts_in.resize(contacts_.size());
  std::copy(contacts_.begin(), contacts_.end(), contacts_in.begin());
}

}  // namespace fcl

#endif
