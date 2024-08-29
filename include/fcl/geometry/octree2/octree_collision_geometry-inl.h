//
// Created by mech-mind_gw on 3/25/2024.
//

#pragma once

#include "fcl/geometry/octree2/octree_prune.h"
#include "fcl/geometry/octree2/octree_visit.h"

namespace fcl {

template <typename S>
Octree2CollisionGeometry<S>::Octree2CollisionGeometry(OctreePtr octree_in)
    : octree(std::move(octree_in)), prune_info(nullptr) {}

template <typename S>
Octree2CollisionGeometry<S>::Octree2CollisionGeometry(
    OctreePtr octree_in, std::shared_ptr<const OctreePruneInfo> prune_info_in)
    : octree(std::move(octree_in)), prune_info(std::move(prune_info_in)) {}

template <typename S>
const std::vector<octree2::OctreeInnerNode>&
Octree2CollisionGeometry<S>::inner_nodes() const {
  return octree->inner_nodes();
}

template <typename S>
const std::vector<bool>&
Octree2CollisionGeometry<S>::inner_nodes_fully_occupied() const {
  if (prune_info != nullptr) {
    return prune_info->new_inner_nodes_fully_occupied;
  }
  return octree->inner_nodes_fully_occupied();
}

template <typename S>
const std::vector<octree2::OctreeLeafNode>&
Octree2CollisionGeometry<S>::leaf_nodes() const {
  if (prune_info != nullptr) {
    return prune_info->new_leaf_nodes;
  }
  return octree->leaf_nodes();
}

template <typename S>
const std::vector<bool>* Octree2CollisionGeometry<S>::prune_internal_nodes()
    const {
  if (prune_info == nullptr) {
    return nullptr;
  }
  return &prune_info->prune_internal_nodes;
}

template <typename S>
const AABB<S>& Octree2CollisionGeometry<S>::octree_root_bv() const {
  return octree->root_bv();
}

template <typename S>
octree2::OctreeTraverseStackElement<S>
Octree2CollisionGeometry<S>::makeStackElementChild(
    const octree2::OctreeTraverseStackElement<S>& parent,
    const AABB<S>& child_AABB, std::uint32_t child_vector_index) const {
  return octree->makeStackElementChild(parent, child_AABB, child_vector_index);
}

template <typename S>
void Octree2CollisionGeometry<S>::computeLocalAABB() {
  if (octree == nullptr) {
    return;
  }

  // Valid tree
  if (octree->n_leaf_nodes() == 0) {
    this->aabb_local = octree->root_bv();
  } else {
    this->aabb_local = octree->leaf_points_AABB();
  }

  // Update center and radius
  this->aabb_center = this->aabb_local.center();
  this->aabb_radius = (this->aabb_local.min_ - this->aabb_center).norm();
}

template <typename S>
void Octree2CollisionGeometry<S>::visitLeafNodes(
    const VisitLeafNodeFunctor& visit_functor) const {
  auto inner_visitor = [&visit_functor](const AABB<S>& aabb, std::uint8_t,
                                        bool is_leaf) -> bool {
    if (!is_leaf) return false;
    return visit_functor(aabb);
  };

  octree2::visitOctree<S>(*octree, prune_info.get(), inner_visitor);
}

template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>>
Octree2CollisionGeometry<S>::pruneBy(const OBB<S>& obb,
                                     bool rebuild_octree) const {
  // Check validity
  if (octree == nullptr) {
    return nullptr;
  }

  // Gather existing prune info
  auto new_prune_info = std::make_shared<OctreePruneInfo>();
  if (prune_info != nullptr) {
    *new_prune_info = *prune_info;
  }

  // Run prune
  octree2::pruneOctreeByOBB(*octree, obb, *new_prune_info);

  // Directly return if do not need rebuild
  if (!rebuild_octree) {
    auto new_tree = std::make_shared<Octree2CollisionGeometry<S>>(
        octree, std::move(new_prune_info));
    new_tree->computeLocalAABB();
    return new_tree;
  }

  // Else, rebuild the tree
  auto new_tree = std::make_shared<octree2::Octree<S>>(*octree);
  new_tree->rebuildAccordingToPruneInfo(*new_prune_info);
  auto new_tree_geom =
      std::make_shared<Octree2CollisionGeometry<S>>(std::move(new_tree));
  new_tree_geom->computeLocalAABB();
  return new_tree_geom;
}

template <typename S>
std::shared_ptr<const Octree2CollisionGeometry<S>>
Octree2CollisionGeometry<S>::rebuildByConsolidatePruneInfo() const {
  auto new_tree = std::make_shared<octree2::Octree<S>>(*octree);
  if (prune_info != nullptr) {
    new_tree->rebuildAccordingToPruneInfo(*prune_info);
  }

  // Into geom
  auto new_tree_geom =
      std::make_shared<Octree2CollisionGeometry<S>>(std::move(new_tree));
  new_tree_geom->computeLocalAABB();
  return new_tree_geom;
}

}  // namespace fcl
