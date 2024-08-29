//
// Created by Wei Gao on 2024/8/16.
//

#pragma once

#include "fcl/broadphase/binary_AABB_tree.h"

namespace fcl {

template <typename S>
using BroadphaseAABB_Tree =
    detail::BinaryAABB_Tree<S, detail::SimpleVectorObjectAllocator>;

}  // namespace fcl
