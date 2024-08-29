//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/common/types.h"

namespace fcl {

enum class TimeOfCollisionRequestType {
  kNotRequested,
  kBoxApproximate,
  kOneTocSample,
  // TODO(wei): support this
  // kExact
};

template <typename S>
struct ContinuousCollisionRequest {
  /// The request for time of collision
  TimeOfCollisionRequestType request_type{
      TimeOfCollisionRequestType::kNotRequested};

  /// The maximum number of contacts that can be returned.
  std::size_t num_max_contacts{1};

  /// Tolerance
  S zero_movement_tolerance{1e-4};
  S gjk_tolerance{1e-6};
  int max_gjk_iterations{128};
};

}  // namespace fcl