//
// Created by mech-mind_gw on 6/18/2022.
//

#pragma once

#include "epa_polytope_utils.h"
#include "minkowski_diff.h"

namespace fcl {
namespace cvx_collide {

enum class GJK_Status {
  // Success case
  // We can certify collision by an simplex that encloses the origin,
  // or certify separation by a direction d such that
  //           dot(d, minkowski_diff.support(d)) < 0
  // If minimum separation is required, a successful estimate is returned
  Intersect,
  Separated,
  // Below as failure/error
  // In these two cases, we are not able to certify collision or separation.
  // An estimate of min-distance is still returned if requested.
  ConvergeNoProgress,
  IterationLimit,
  // Unknown failure
  Failed
};

template <typename T>
struct GJKSimplex {
  MinkowskiDiffVertex<T> vertices[4];
  int rank{-1};

  // Access
  void reset() { rank = -1; }
  bool is_valid() const { return rank > 0; }
  int n_vertices() const { return rank; }

  // Add vertex
  void AddVertex(MinkowskiDiffVertex<T> vertex) {
    if (rank < 0) rank = 0;
    assert(rank < 4);
    vertices[rank] = std::move(vertex);
    rank += 1;
  }
};

template <typename T>
class GJK {
 public:
  GJK(std::size_t max_iterations, T tolerance)
      : max_iterations_(max_iterations), tolerance_(tolerance) {}

  // If separation distance is requested, the output is written into this
  // struct which contains the two point that reach the min separation distance
  // on geometry_0 and geometry_1, BOTH EXPRESSED IN geometry_0 FRAME
  struct MinSeparationDistanceOutput {
    // The vertex such that dot(d, minkowski_diff.support(d)) < 0
    // which certify separation
    bool is_certification_vertex_valid{false};
    MinkowskiDiffVertex<T> vertex_certify_separation;

    // The two separating points that achieve minimum distance
    bool is_separation_point_valid{false};
    Vector3<T> p0_if_separated;
    Vector3<T> p1_if_separated;

    // Init as invalid
    MinSeparationDistanceOutput() = default;

    // clang-format off
    Vector3<T> point_on_minkowskidiff() const { return p0_if_separated - p1_if_separated; }
    const Vector3<T>& p0_on_shape0_frame() const { return p0_if_separated; }
    const Vector3<T>& p1_on_shape0_frame() const { return p1_if_separated; }
    T separation_distance() const { return (p0_if_separated - p1_if_separated).norm(); }
    T separation_distance_squared() const { return (p0_if_separated - p1_if_separated).squaredNorm(); }
    // clang-format on
  };

  // Run interface
  // The simplex as an output depends on the returned status
  // If status == Intersect, then simplex contains the origin
  // Note that there might be touching containment
  // If status == Separated, the simplex is the converging one
  GJK_Status Evaluate(const MinkowskiDiff<T>& shape, GJKSimplex<T>& simplex,
                      const Vector3<T>& guess = Vector3<T>::UnitX(),
                      MinSeparationDistanceOutput*
                          min_distance_output_if_separated = nullptr) const;

 private:
  // Parameters
  const std::size_t max_iterations_;
  const T tolerance_;

  // The projection interface
  enum class ProjectionStatus {
    Failed,
    Continue,
    Intersect,
    FailedTetrahedronZeroVolume
  };
  ProjectionStatus simplexProjection(GJKSimplex<T>& simplex,
                                     Vector3<T>& direction) const;
  ProjectionStatus simplexProjection2(GJKSimplex<T>& simplex,
                                      Vector3<T>& direction) const;
  ProjectionStatus simplexProjection3(GJKSimplex<T>& simplex,
                                      Vector3<T>& direction) const;
  ProjectionStatus simplexProjection4(GJKSimplex<T>& simplex,
                                      Vector3<T>& direction) const;

  // The separation distance interface
  // simplex is both input and output, as input it contains ONE vertex
  // that witness the separation of the point
  // return:
  bool findMinimumDistancePointsWithSeparatedVertexInit(
      const MinkowskiDiff<T>& shape, GJKSimplex<T>& simplex,
      std::pair<Vector3<T>, Vector3<T>>& p0p1_in_frame0) const;

  // The functions for update the min-distance simplex, although
  // the last vertex in the simplex is the latest one, this information
  // can only be used as speedup.
  // If return status == MinDistanceUpdateStatus::OK:
  //    simplex are updated, min_distance_output is meaningful
  //    However, it might be no improvement that we CANNOT detect.
  // If return status == NoImprovement:
  //    simplex are updated to original (without new vertex)
  //    min_distance_output is NOT meaningful
  //    Direct return after this method
  enum MinDistanceUpdateStatus { NoImprovement, OK, Failed };
  static constexpr T barycentric_weight_tolerance = 1e-3;
  MinDistanceUpdateStatus computeMinDistanceAndUpdateSimplex(
      GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const;
  // There status can ONLY be OK in simplex2/3
  void computeMinDistanceAndUpdateSimplex2(
      GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const;
  void computeMinDistanceAndUpdateSimplex3(
      GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const;
  MinDistanceUpdateStatus computeMinDistanceAndUpdateSimplex4(
      GJKSimplex<T>& simplex, Vector3<T>& min_distance_output) const;

  // Extract the separation point in the given simplex
  // Return: whether the min-distance point is on the given simplex
  //         or return output.is_separation_point_valid
  // If return false, elements in the output are usually not valid.
  // However, to ease the implementation when simplex.size == 2
  // (segment simplex) the min distance point is always write into the
  // output (even if the min distance is not achieved on segment, but on vertex)
  bool extractSeparationPointNoSubSimplex(
      const MinkowskiDiff<T>& shape, const GJKSimplex<T>& simplex,
      std::pair<Vector3<T>, Vector3<T>>& p0p1_in_frame0) const;
  bool extractSeparationPointTrySubSimplex(
      const MinkowskiDiff<T>& shape, const GJKSimplex<T>& simplex,
      std::pair<Vector3<T>, Vector3<T>>& p0p1_in_frame0) const;
};

}  // namespace cvx_collide
}  // namespace fcl

#include "gjk.hpp"
#include "gjk_distance.hpp"
