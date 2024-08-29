//
// Created by mech-mind_gw on 12/25/2023.
//

#pragma once

#include <cassert>
#include <unordered_set>

#include "fcl/geometry/shape/triangle_p.h"

namespace fcl {

template <typename S>
void makeTriangleSoupSurfaceVoxel(
    const std::vector<TriangleP<S>>& triangles, S resolution,
    std::unordered_set<std::int64_t>& encoded_voxel_set);

}  // namespace fcl

#include "octree_from_triangles-inl.h"