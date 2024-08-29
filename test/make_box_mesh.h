#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/math/mesh_simplex.h"

namespace fcl {
namespace test {

/*
 Generates unique vertices on a Cartesian grid of the box. In each of the
 x-, y-, and z-directions, the vertices are distributed uniformly.
 @param[in] box
     The box shape specification (see drake::geometry::Box).
 @param[in] num_vertices
     Number of vertices in each of x-, y-, and z-directions.
 @retval vertices
     The linear sequence of vertices consistent with CalcSequentialIndex.
 */
template <typename T>
std::vector<Vector3<T>> GenerateVertices(const Box<T>& box,
                                         const Vector3<int>& num_vertices) {
  const T half_x = box.side[0] / T(2);
  const T half_y = box.side[1] / T(2);
  const T half_z = box.side[2] / T(2);
  const auto x_coords =
      VectorX<T>::LinSpaced(num_vertices.x(), -half_x, half_x);
  const auto y_coords =
      VectorX<T>::LinSpaced(num_vertices.y(), -half_y, half_y);
  const auto z_coords =
      VectorX<T>::LinSpaced(num_vertices.z(), -half_z, half_z);

  std::vector<Vector3<T>> vertices;
  vertices.reserve(num_vertices.x() * num_vertices.y() * num_vertices.z());
  // The order of nested i-loop, j-loop, then k-loop makes the sequence of
  // vertices consistent with CalcSequentialIndex.
  for (int i = 0; i < num_vertices.x(); ++i) {
    for (int j = 0; j < num_vertices.y(); ++j) {
      for (int k = 0; k < num_vertices.z(); ++k) {
        vertices.emplace_back(x_coords[i], y_coords[j], z_coords[k]);
      }
    }
  }
  return vertices;
}

/*
 Calculates the sequential vertex index of the vertex specified by the
 (i, j, k)-index in the Cartesian grid defined by the number of vertices in
 x-, y-, and z-directions.
 @param i  Index in x-direction.
 @param j  Index in y-direction.
 @param k  Index in z-direction.
 @param num_vertices  Number of vertices in x-, y-, and z-directions.
 @pre 0 ≤ i < num_vertices.x(),
      0 ≤ j < num_vertices.y(), and
      0 ≤ k < num_vertices.z().
 */
int CalcSequentialIndex(int i, int j, int k, const Vector3<int>& num_vertices) {
  return i * num_vertices.y() * num_vertices.z() + j * num_vertices.z() + k;
}

/*
 Adds six tetrahedra of a given rectangular cell to the list of tetrahedral
 elements. The rectangular cell is identified by the (i,j,k)-index of its
 lowest vertex. The (i,j,k)-indices of its 8 vertices are
 (i,j,k) + {0,1}x{0,1}x{0,1}.
 @param[in] lowest
     The (i,j,k) index of the lowest vertex of the rectangular cell.
 @param[in] num_vertices
     Number of vertices in each of x-, y-, and z-directions.
 @param[in,out] elements
     The six tetrahedra are added into this list of elements.
 */
void AddSixTetrahedraOfCell(const Vector3<int>& lowest,
                            const Vector3<int>& num_vertices,
                            std::vector<MeshSimplex>* elements) {
  const int i = lowest.x();
  const int j = lowest.y();
  const int k = lowest.z();
  // Get the sequential indices of the eight vertices of the rectangular cell.
  int v[8];
  int s = 0;
  for (int l = 0; l < 2; ++l)
    for (int m = 0; m < 2; ++m)
      for (int n = 0; n < 2; ++n)
        v[s++] = CalcSequentialIndex(i + l, j + m, k + n, num_vertices);
  // The following picture shows where vertex vₛ (for `v[s]` above) locates
  // in the rectangular cell.  The I-, J-, K-axes show the direction of
  // increasing i, j, k indices.
  //
  //               v₁     v₃
  //               ●------●
  //              /|     /|
  //             / |  v₇/ |
  //         v₅ ●------●  |
  //            |  |   |  |
  //            |  ●---|--● v₂
  //            | /v₀  | /
  //            |/     |/
  //    +K   v₄ ●------● v₆
  //     |
  //     |
  //     o------+J
  //    /
  //   /
  // +I
  //
  // The following table subdivides the rectangular cell into six tetrahedra.
  // Refer to the picture above to determine which four vertices form a
  // tetrahedron. The six tetrahedra form a cycle around the diagonal v₀v₇
  // of the cell. Refer to:
  // http://www.baumanneduard.ch/Splitting%20a%20cube%20in%20tetrahedras2.htm
  const int tetrahedron[6][4]{
      // clang-format off
      {v[0], v[7], v[4], v[6]},
      {v[0], v[7], v[6], v[2]},
      {v[0], v[7], v[2], v[3]},
      {v[0], v[7], v[3], v[1]},
      {v[0], v[7], v[1], v[5]},
      {v[0], v[7], v[5], v[4]}};
  // clang-format on

  // The above table guarantees that adjacent rectangular cells will be
  // subdivided in a consistent way, i.e., both cells will pick the same
  // diagonal of their shared rectangular face.
  for (int t = 0; t < 6; ++t)
    elements->emplace_back(MeshSimplex(tetrahedron[t][0], tetrahedron[t][1],
                                       tetrahedron[t][2], tetrahedron[t][3]));
}

/*
 Generates connectivity for the tetrahedral elements of the mesh.
 @param[in] num_vertices
     Number of vertices in each of x-, y-, and z-directions.
 @return
     A sequence of tetrahedral elements that share unique vertices.
 */
std::vector<MeshSimplex> GenerateElements(const Vector3<int>& num_vertices) {
  std::vector<MeshSimplex> elements;
  const Vector3<int> num_cell = num_vertices - Vector3<int>(1, 1, 1);
  elements.reserve(6 * num_cell.x() * num_cell.y() * num_cell.z());
  for (int i = 0; i < num_cell.x(); ++i) {
    for (int j = 0; j < num_cell.y(); ++j) {
      for (int k = 0; k < num_cell.z(); ++k) {
        AddSixTetrahedraOfCell(Vector3<int>(i, j, k), num_vertices, &elements);
      }
    }
  }
  return elements;
}

/*
 Generates a tetrahedral volume mesh of a given box by subdividing the box
 into _rectangular cells_ (volume bounded by six axis-aligned faces) and
 subdividing each rectangular cell into six tetrahedra. The output mesh will
 have these properties:

   1. The generated vertices are unique. There is no repeating vertices in
      the list of vertex coordinates.
   2. The generated tetrahedra are _conforming_. Two tetrahedra intersect in
      their shared face, or shared edge, or shared vertex, or not at all.
      There is no partial overlapping of two tetrahedra.

 @param[in] box
     The box shape specification (see drake::geometry::Box).
 @param[in] resolution_hint
     A measure (in meters) that controls the resolution of the mesh. The length
     of the axis-aligned edges of the mesh will be less than or equal to this
     parameter. The length of non-axis-aligned edges will be less than or equal
     to √3 of this parameter. The coarsest possible mesh can be made by
     providing a resolution hint at least as large as the box's largest
     dimension.
 @retval volume_mesh
 @tparam_nonsymbolic_scalar
 @note The mesh has no guarantee on the inner boundary for a rigid core.
 */
template <typename BV>
BVHModel<BV> MakeBoxBVHTetrahedronModel(const Box<typename BV::S>& box,
                                        double resolution_hint) {
  using S = typename BV::S;
  // Number of vertices in x-, y-, and z- directions.  In each direction,
  // there is one more vertices than cells.
  const Vector3<int> num_vertices{
      1 + static_cast<int>(ceil(box.side[0] / resolution_hint)),
      1 + static_cast<int>(ceil(box.side[1] / resolution_hint)),
      1 + static_cast<int>(ceil(box.side[2] / resolution_hint))};

  const std::vector<Vector3<S>> vertices =
      GenerateVertices<S>(box, num_vertices);

  const std::vector<MeshSimplex> elements = GenerateElements(num_vertices);

  BVHModel<BV> bvh;
  bvh.beginModel();
  bvh.addSubModel(vertices, elements);
  bvh.endModel();
  return bvh;
}

}  // namespace test
}  // namespace fcl
