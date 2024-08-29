//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

#include <functional>
#include <stack>
#include <vector>

#include "fcl/geometry/octree2/octree.h"

namespace fcl {
namespace octree2 {

/// Visitor function. The is_leaf parameter here implies whether the currently
/// visited node has any children. Do NOT confused it with the LeafNode or
/// RealNodeNoLeaf that is INTERNALLY defined for octree.
template <typename S>
using VisitOctreeNodeFunc =
    std::function<bool(const AABB<S>& aabb, std::uint8_t depth, bool is_leaf)>;

/// Actual visiting routines
template <typename S>
void visitOctree(const Octree<S>& tree,
                 const OctreePruneInfo* prune_octree_info,
                 const VisitOctreeNodeFunc<S>& visitor);
template <typename S>
void visitOctree(const Octree<S>& tree, const VisitOctreeNodeFunc<S>& visitor);

}  // namespace octree2
}  // namespace fcl

#include "fcl/geometry/octree2/octree_visit-inl.h"
