//
// Created by mech-mind_gw on 12/16/2023.
//

#pragma once

namespace fcl {
namespace detail {

template <typename S>
inline void find_min_max(S x0, S x1, S x2, S& min, S& max);
template <typename S>
inline void find_min_max(S x0, S x1, S x2, S x3, S& min, S& max);

}
}  // namespace fcl

#include "fcl/narrowphase/detail/primitive_shape_algorithm/simplex_common-inl.h"