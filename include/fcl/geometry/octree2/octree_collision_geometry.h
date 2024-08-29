//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

#include "fcl/geometry/collision_geometry.h"
#include "fcl/geometry/octree2/octree.h"

namespace fcl {

template <typename S>
class Octree2CollisionGeometry : public CollisionGeometry<S> {
 public:
  using Octree = octree2::Octree<S>;
  using OctreePtr = std::shared_ptr<const Octree>;
  using OctreeInnerNode = octree2::OctreeInnerNode;
  using OctreeLeafNode = octree2::OctreeLeafNode;
  using OctreeTraverseStackElement = octree2::OctreeTraverseStackElement<S>;
  using OctreePruneInfo = octree2::OctreePruneInfo;
  using PruneInfoPtr = std::shared_ptr<const OctreePruneInfo>;
  using ConstPtr = std::shared_ptr<const Octree2CollisionGeometry<S>>;
  explicit Octree2CollisionGeometry(OctreePtr octree);
  Octree2CollisionGeometry(OctreePtr octree,
                           std::shared_ptr<const OctreePruneInfo> prune_info);
  ~Octree2CollisionGeometry() override = default;

  /// Compute the AABB<S> for the geometry in its local coordinate system
  void computeLocalAABB() override;

  /// Remove an OBB from this tree
  ConstPtr pruneBy(const OBB<S>& obb, bool rebuild_octree) const;
  ConstPtr rebuildByConsolidatePruneInfo() const;

  /// Simple access
  // clang-format off
  OBJECT_TYPE getObjectType() const override { return OT_OCTREE2; };
  NODE_TYPE getNodeType() const override { return GEOM_OCTREE2; };
  const OctreePtr& raw_octree() const { return octree; };
  const std::vector<OctreeInnerNode>& inner_nodes() const;
  const std::vector<bool>& inner_nodes_fully_occupied() const;
  const std::vector<OctreeLeafNode>& leaf_nodes() const;
  const std::vector<bool>* prune_internal_nodes() const;
  const AABB<S>& octree_root_bv() const;
  OctreeTraverseStackElement makeStackElementChild(
    const OctreeTraverseStackElement& parent, const AABB<S>& child_AABB,
    std::uint32_t child_vector_index) const;
  // clang-format on

  /// Visit all the bbox of the leaf nodes with a user-provided
  /// functor. c/t are the uncertainty information.
  using VisitLeafNodeFunctor = std::function<bool(const AABB<S>& aabb)>;
  void visitLeafNodes(const VisitLeafNodeFunctor& visit_functor) const;

  /// Read-only, shared ptr access to a octree with pruned or not
 private:
  std::shared_ptr<const Octree> octree{nullptr};
  std::shared_ptr<const OctreePruneInfo> prune_info{nullptr};
};

}  // namespace fcl

#include "fcl/geometry/octree2/octree_collision_geometry-inl.h"
