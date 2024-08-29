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

#ifndef FCL_COLLISIONRESULT_H
#define FCL_COLLISIONRESULT_H

#include <set>
#include <vector>

#include "fcl/common/types.h"
#include "fcl/narrowphase/contact.h"

namespace fcl {

/// @brief user can define a functor the process the collision contact
///        and determine:
///        1. Shall we keep this contact
///        2. Can we terminate the narrowphase collision
///        Note that even through the user suggest termination, the algorithm
///        may not terminate immediately
template <typename S>
using UserContactProcessFunctor =
    std::function<void(const Contact<S>& c, bool& output_keep_this_contact,
                       bool& output_can_we_terminate)>;

/// @brief collision result
template <typename S>
struct CollisionResult {
 private:
  /// @brief contact information
  std::vector<Contact<S>> contacts_;

  /// @brief optional user specified functor
  UserContactProcessFunctor<S> user_process_functor_;
  bool user_stop_{false};

 public:
  explicit CollisionResult();
  explicit CollisionResult(UserContactProcessFunctor<S> user_functor);
  ~CollisionResult() = default;

  /// @brief The interface implementation
  void addContact(const Contact<S>& c);
  std::size_t numContacts() const;
  void clear();

  /// @brief Binary collision checking, by default it is as numContacts() > 0
  bool isCollision() const;

  /// @brief Determine whether the collision checking can stop
  ///        Default implementation match the behavior of fcl
  bool terminationConditionSatisfied(std::size_t n_max_contacts) const;

  /// @brief Get the explicitly stored contacts
  const Contact<S>& getContact(std::size_t i) const;
  const std::vector<Contact<S>>& getContacts() const { return contacts_; }
  std::vector<Contact<S>>& getContacts() { return contacts_; }
  void getContacts(std::vector<Contact<S>>& contacts);
};

using CollisionResultf = CollisionResult<float>;
using CollisionResultd = CollisionResult<double>;

}  // namespace fcl

#include "fcl/narrowphase/collision_result-inl.h"

#endif
