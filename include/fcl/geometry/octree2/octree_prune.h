//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_util.h"

namespace fcl {
namespace octree2 {

/// Remove an obb from the Octree, the tree might be already pruned by this
/// method. If so, existing_prune_info might be valid.
template <typename S>
void pruneOctreeByOBB(const Octree<S>& tree, const OBB<S>& pruned_obb,
                      OctreePruneInfo& new_or_existing_prune_info);

}  // namespace octree2
}  // namespace fcl

#include "fcl/geometry/octree2/octree_prune-inl.h"
