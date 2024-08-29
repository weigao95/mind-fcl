//
// Created by Wei Gao on 2024/8/15.
//
#include "fcl/broadphase/binary_AABB_tree.h"

namespace fcl {
namespace detail {

template class BinaryAABB_Tree<float, SimpleVectorObjectAllocator>;
template class BinaryAABB_Tree<double, SimpleVectorObjectAllocator>;

}  // namespace detail
}  // namespace fcl
