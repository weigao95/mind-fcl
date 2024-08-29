//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

#include "fcl/geometry/octree2/octree_node.h"
#include "fcl/math/bv/AABB.h"
#include "fcl/math/bv/OBB.h"

namespace fcl {
namespace octree2 {

/// Compute the bounding volume of an octree node's i-th child
template <typename S>
void computeChildAABB(const AABB<S>& parent_bv, std::uint8_t child_i,
                      AABB<S>& child_bv);

/// Traversing stack element for visiting a octree
template <typename S>
struct OctreeTraverseStackElement {
  AABB<S> bv;
  std::uint8_t depth;
  bool is_leaf_node;
  std::uint32_t node_vector_index;

  static OctreeTraverseStackElement<S> MakeRoot(const AABB<S>& root_bv);
};

/// Containment test of AABB and OBB
template <typename S>
bool is_contained_naive(const OBB<S>& container,
                        const AABB<S>& maybe_contained);
template <typename S>
bool is_contained(const OBB<S>& container, const AABB<S>& maybe_contained);

}  // namespace octree2
}  // namespace fcl

#include "fcl/geometry/octree2/octree_util-inl.h"
