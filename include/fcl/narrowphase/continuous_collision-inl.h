#pragma once

#include "fcl/narrowphase/detail/ccd/translational_collision_func_matrix.h"

namespace fcl {
namespace detail {

//==============================================================================
template <typename S>
const TranslationalCollisionFunctionMatrix<S>&
translationalCollisionFunctionMatrix() {
  static TranslationalCollisionFunctionMatrix<S> table;
  return table;
}

}  // namespace detail
}  // namespace fcl

namespace fcl {

template <typename S>
void translational_ccd(const CollisionGeometry<S>* o1, const Transform3<S>& tf1,
                       const TranslationalDisplacement<S>& o1_displacement,
                       const CollisionGeometry<S>* o2, const Transform3<S>& tf2,
                       const ContinuousCollisionRequest<S>& request,
                       ContinuousCollisionResult<S>& result) {
  const auto& look_table = detail::translationalCollisionFunctionMatrix<S>();
  const auto node_type1 = o1->getNodeType();
  const auto node_type2 = o2->getNodeType();
  if (!look_table.collision_matrix[node_type1][node_type2]) {
    std::cerr << "Warning: collision function between node type " << node_type1
              << " and node type " << node_type2 << " is not supported\n";
  } else {
    look_table.collision_matrix[node_type1][node_type2](
        o1, tf1, o1_displacement, o2, tf2, request, result);
  }
}

}  // namespace fcl