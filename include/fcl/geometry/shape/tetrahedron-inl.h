#pragma once

namespace fcl {

template <typename S>
Tetrahedron<S>::Tetrahedron(const Vector3<S>& v0, const Vector3<S>& v1,
                            const Vector3<S>& v2, const Vector3<S>& v3)
    : Tetrahedron(
          std::array<Vector3<S>, kVertexOfTetrahedron>{v0, v1, v2, v3}) {}

/// @brief virtual function of compute AABB<S> in local coordinate
template <typename S>
void Tetrahedron<S>::computeLocalAABB() {
  this->aabb_local.min_.setConstant(std::numeric_limits<S>::max());
  this->aabb_local.max_.setConstant(-std::numeric_limits<S>::max());

  for (const auto& v : vertices) {
    this->aabb_local += v;
  }

  this->aabb_center = this->aabb_local.center();
  this->aabb_radius = (this->aabb_local.min_ - this->aabb_center).norm();
}

// Documentation inherited
template <typename S>
NODE_TYPE Tetrahedron<S>::getNodeType() const {
  return GEOM_TETRAHEDRON;
}

template <typename S>
std::vector<Vector3<S>> Tetrahedron<S>::getBoundVertices(
    const Transform3<S>& tf) const {
  std::vector<Vector3<S>> result;
  result.reserve(4);
  for (const auto& vert : vertices) {
    result.emplace_back(tf * vert);
  }
  return result;
}

}  // namespace fcl