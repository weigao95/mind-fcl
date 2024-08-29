//
// Created by Wei Gao on 2024/8/16.
//

#pragma once

#include "fcl/common/types.h"
#include "fcl/math/bv/AABB.h"

namespace fcl {

template <typename S>
struct BroadphaseObjectInfo {
  AABB<S> bv{};
  std::uint64_t user_id{0};

  // Constructors
  BroadphaseObjectInfo() = default;
  BroadphaseObjectInfo(AABB<S> aabb, std::uint64_t user_id_in)
      : bv(std::move(aabb)), user_id(user_id_in) {}
};

}  // namespace fcl
