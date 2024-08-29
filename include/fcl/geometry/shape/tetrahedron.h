//
// Created by mech-mind_gw on 10/21/2023.
//

#pragma once

#include "fcl/geometry/shape/shape_base.h"

namespace fcl {


/// @brief Tetrahedron with 4 points. This class would be used internally (in
/// the same way as TriangleP). There won't be a fcl::collide interface for it.
template <typename S>
class Tetrahedron : public ShapeBase<S> {
 public:
  /// @brief Default constructor
  static constexpr int kVertexOfTetrahedron = 4;
  Tetrahedron() = default;

  explicit Tetrahedron(std::array<Vector3<S>, 4> vertices_in)
      : vertices(std::move(vertices_in)) {}

  /// @brief Constructor with individual vertices
  /// @param v0,v1,v2,v3 Vertices of the tetrahedron
  Tetrahedron(const Vector3<S>& v0, const Vector3<S>& v1, const Vector3<S>& v2,
              const Vector3<S>& v3);

  /// @brief virtual function of compute AABB<S> in local coordinate
  void computeLocalAABB() override;

  // Documentation inherited
  NODE_TYPE getNodeType() const override;

  /// @brief Array to store each vertex of the tetrahedron
  std::array<Vector3<S>, kVertexOfTetrahedron> vertices;

  /// @brief get the vertices of some convex shape which can bound this shape in
  /// a specific configuration
  std::vector<Vector3<S>> getBoundVertices(const Transform3<S>& tf) const;
};

}  // namespace fcl

#include "fcl/geometry/shape/tetrahedron-inl.h"