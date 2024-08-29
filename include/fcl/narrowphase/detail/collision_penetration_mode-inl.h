//
// Created by mech-mind_gw on 3/31/2023.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
CollisionPenetrationMode<S>::CollisionPenetrationMode()
    : type_(CollisionPenetrationType::Disabled) {
  unit_direction_data_ = Vector3<S>::UnitZ();
}

template <typename S>
CollisionPenetrationMode<S>::CollisionPenetrationMode(
    CollisionPenetrationType type)
    : type_(type) {
  unit_direction_data_ = Vector3<S>::UnitZ();
}

template <typename S>
CollisionPenetrationMode<S>::CollisionPenetrationMode(
    CollisionPenetrationType type, Vector3<S> unit_direction)
    : type_(type), unit_direction_data_(std::move(unit_direction)) {
  assert(std::abs(unit_direction_data_.norm() - S(1.0)) < S(1e-3));
}

template <typename S>
void CollisionPenetrationMode<S>::setAsDisabled() {
  type_ = CollisionPenetrationType::Disabled;
}

template <typename S>
void CollisionPenetrationMode<S>::setAsDefaultGJK_EPA() {
  type_ = CollisionPenetrationType::DefaultGJK_EPA;
}

template <typename S>
void CollisionPenetrationMode<S>::setAsDirectedPenetration(
    Vector3<S> request_direction_unit) {
  type_ = CollisionPenetrationType::DirectedPenetration;
  unit_direction_data_ = std::move(request_direction_unit);
  assert(std::abs(unit_direction_data_.norm() - S(1.0)) < S(1e-3));
}

template <typename S>
void CollisionPenetrationMode<S>::setAsIncrementalMinimumPenetration(
    Vector3<S> direction_hint_unit) {
  type_ = CollisionPenetrationType::IncrementalMinimumPenetration;
  unit_direction_data_ = std::move(direction_hint_unit);
  assert(std::abs(unit_direction_data_.norm() - S(1.0)) < S(1e-3));
}

template <typename S>
const Vector3<S>& CollisionPenetrationMode<S>::priorPenetrationDirection()
    const {
  assert(type_ == CollisionPenetrationType::DirectedPenetration ||
         type_ == CollisionPenetrationType::IncrementalMinimumPenetration);
  return unit_direction_data_;
}

// common types
template <typename S>
const CollisionPenetrationMode<S>& CollisionPenetrationMode<S>::Disabled() {
  static CollisionPenetrationMode<S> disabled(
      CollisionPenetrationType::Disabled);
  return disabled;
}

template <typename S>
const CollisionPenetrationMode<S>&
CollisionPenetrationMode<S>::MinimumPenetrationEPA() {
  static CollisionPenetrationMode<S> disabled(
      CollisionPenetrationType::DefaultGJK_EPA);
  return disabled;
}

}  // namespace detail
}  // namespace fcl
