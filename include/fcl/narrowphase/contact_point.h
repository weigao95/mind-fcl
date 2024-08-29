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

#ifndef FCL_CONTACTPOINT_H
#define FCL_CONTACTPOINT_H

#include "fcl/common/types.h"

namespace fcl {

/// @brief Minimal contact information returned by collision. Only holds the
///        penetration contact (i.e., no separation distance).
///
///        Another method to represent the contact is by two points p1/p2 on
///        geometry 1/2. These two points constitute the penetration distance.
///        Note that when two objects overlap, the direction
///                  direction_point_1_to_2 = p2 - p1
///        can be OPPOSITE to the intuitive "1_to_2" direction, which is
///        usually defined as
///                  direction_intuitive_1_to_2 = center_2 - center_1
///        where center_1/2 are the geometric center of corresponded geometry.
///
///        Besides, there is also an "escape" direction in motion planning.
///        If geom2 is moved in
///                  direction_escape_1_to_2 = p1 - p2
///        for distance equals (positive) penetration_depth, the geom2 escapes
///        from colliding with geom1. In other words, the direction and
///        distance such that point p2 on geom2 move to escape from geom1.
///        Using this definition, the direction we would use (normal) in the
///        class is "direction_escape_1_to_2", this is roughly in the same
///        direction as "direction_intuitive_1_to_2", and is the NEGATIVE of
///        "direction_point_1_to_2".
/// \tparam S float/double
template <typename S>
struct ContactPoint {
  /// @brief Contact normal, the direction_escape_1_to_2 mentioned above.
  ///        normalized(p1 - p2);
  Vector3<S> normal;

  /// @brief Contact position, in world space 0.5 * (p1 + p2)
  Vector3<S> pos;

  /// @brief Penetration depth, always POSITIVE despite this definition is
  ///        not aligned with the convention signed distance.
  S penetration_depth;

  /// @brief Constructor
  ContactPoint();

  /// @brief Constructor
  ContactPoint(const Vector3<S>& normal_in, const Vector3<S>& position_in,
               S penetration_depth_in);

  /// Obtain the original point on p1 and p2
  Vector3<S> point_on_shape1() const;
  Vector3<S> point_on_shape2() const;
  Vector3<S> shape2_escape_movement() const;
};

/// @brief Array-based contact information returned by collision
template <typename S>
class CollisionPenetrationContactData {
 private:
  /// In narrowphase, the max number of contact point is box2Box
  /// There can be at most 8, but it is restricted to 4 in our setups
  static constexpr std::size_t kNarrowphaseMaxNumContacts = 4;

  // Actual container and size index
  using ArrayType = std::array<ContactPoint<S>, kNarrowphaseMaxNumContacts>;
  ArrayType contacts_;
  std::size_t n_contacts_{0};

 public:
  // Default construct
  explicit CollisionPenetrationContactData() = default;
  ~CollisionPenetrationContactData() = default;

  // Default copy/assign/move
  // clang-format off
  CollisionPenetrationContactData(const CollisionPenetrationContactData<S>&) = default;
  CollisionPenetrationContactData(CollisionPenetrationContactData<S>&&) noexcept = default;
  CollisionPenetrationContactData& operator=(const CollisionPenetrationContactData<S>&) = default;
  CollisionPenetrationContactData& operator=(CollisionPenetrationContactData<S>&&) noexcept = default;
  // clang-format on

  // Simple query
  std::size_t size() const { return n_contacts_; }
  bool empty() const { return n_contacts_ == 0; }
  void clear() { n_contacts_ = 0; }
  void resize(std::size_t new_size) {
    assert(new_size <= kNarrowphaseMaxNumContacts);
    n_contacts_ = new_size;
  }

  // Operators
  ContactPoint<S>& operator[](std::size_t i) { return contacts_[i]; }
  const ContactPoint<S>& operator[](std::size_t i) const {
    return contacts_[i];
  }

  // Iterators
  typename ArrayType::iterator begin() { return contacts_.begin(); }
  typename ArrayType::iterator end() { return contacts_.begin() + n_contacts_; }
  typename ArrayType::const_iterator begin() const { return contacts_.begin(); }
  typename ArrayType::const_iterator end() const {
    return contacts_.begin() + n_contacts_;
  }

  // Emplace back
  template <class... Args>
  void emplace_back(Args&&... args) {
    assert(n_contacts_ < contacts_.size());
    if (n_contacts_ >= contacts_.size()) return;
    new (contacts_.data() + n_contacts_)
        ContactPoint<S>(std::forward<Args>(args)...);
    n_contacts_ += 1;
  }
};

/// @brief Return true if _cp1's penetration depth is less than _cp2's.
template <typename S>
bool comparePenDepth(const ContactPoint<S>& _cp1, const ContactPoint<S>& _cp2);

template <typename S>
void flipNormal(CollisionPenetrationContactData<S>& contacts);

}  // namespace fcl

#include "fcl/narrowphase/contact_point-inl.h"

#endif
