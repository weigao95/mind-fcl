//
// Created by Wei Gao on 2024/3/24.
//

#pragma once

namespace fcl {
namespace octree2 {

template <typename S>
void computeChildAABB(const AABB<S>& root_bv, std::uint8_t child_i,
                      AABB<S>& child_bv) {
  constexpr std::uint8_t bit_1 = 1;
  constexpr std::uint8_t bit_10 = 2;
  constexpr std::uint8_t bit_100 = 4;
  if (child_i & bit_1) {
    child_bv.min_[0] = (root_bv.min_[0] + root_bv.max_[0]) * 0.5;
    child_bv.max_[0] = root_bv.max_[0];
  } else {
    child_bv.min_[0] = root_bv.min_[0];
    child_bv.max_[0] = (root_bv.min_[0] + root_bv.max_[0]) * 0.5;
  }

  if (child_i & bit_10) {
    child_bv.min_[1] = (root_bv.min_[1] + root_bv.max_[1]) * 0.5;
    child_bv.max_[1] = root_bv.max_[1];
  } else {
    child_bv.min_[1] = root_bv.min_[1];
    child_bv.max_[1] = (root_bv.min_[1] + root_bv.max_[1]) * 0.5;
  }

  if (child_i & bit_100) {
    child_bv.min_[2] = (root_bv.min_[2] + root_bv.max_[2]) * 0.5;
    child_bv.max_[2] = root_bv.max_[2];
  } else {
    child_bv.min_[2] = root_bv.min_[2];
    child_bv.max_[2] = (root_bv.min_[2] + root_bv.max_[2]) * 0.5;
  }
}

template <typename S>
OctreeTraverseStackElement<S> OctreeTraverseStackElement<S>::MakeRoot(
    const AABB<S>& root_bv) {
  OctreeTraverseStackElement root;
  root.bv = root_bv;
  root.depth = 0;
  root.is_leaf_node = false;
  root.node_vector_index = 0;
  return root;
}

template <typename S>
bool is_contained_naive(const OBB<S>& container,
                        const AABB<S>& maybe_contained) {
  const Vector3<S> center = maybe_contained.center();
  const Vector3<S> offset =
      S(0.5) * (maybe_contained.max_ - maybe_contained.min_);
  const std::array<S, 2> positive_and_negative{-1.0, 1.0};
  Vector3<S> signed_offset;
  for (const auto x : positive_and_negative) {
    for (const auto y : positive_and_negative) {
      for (const auto z : positive_and_negative) {
        signed_offset = offset;
        signed_offset.x() *= x;
        signed_offset.y() *= y;
        signed_offset.z() *= z;
        const Vector3<S> point = center + signed_offset;
        if (!container.contain(point)) {
          return false;
        }
      }
    }
  }

  // All points are tested
  return true;
}

template <typename S>
bool is_contained(const OBB<S>& container, const AABB<S>& maybe_contained) {
  const Vector3<S> aabb_center = maybe_contained.center();
  const Vector3<S> offset =
      S(0.5) * (maybe_contained.max_ - maybe_contained.min_);
  const Vector3<S> aabb_center_in_obb_frame =
      container.axis.transpose() * (aabb_center - container.To);

  // For each obb axis
  for(auto obb_axis_j = 0; obb_axis_j < 3; obb_axis_j++) {
    S axis_min = aabb_center_in_obb_frame[obb_axis_j];
    S axis_max = aabb_center_in_obb_frame[obb_axis_j];
    const auto obb_axis = container.axis.col(obb_axis_j);
    for(auto aabb_axis_i = 0; aabb_axis_i < 3; aabb_axis_i++) {
      if (obb_axis[aabb_axis_i] > 0) {
        axis_max += obb_axis[aabb_axis_i] * offset[aabb_axis_i];
        axis_min -= obb_axis[aabb_axis_i] * offset[aabb_axis_i];
      } else {
        axis_max -= obb_axis[aabb_axis_i] * offset[aabb_axis_i];
        axis_min += obb_axis[aabb_axis_i] * offset[aabb_axis_i];
      }
    }

    if (axis_max > container.extent[obb_axis_j]) return false;
    if (axis_min < -container.extent[obb_axis_j]) return false;
  }

  // Done
  return true;
}

}  // namespace octree2
}  // namespace fcl
