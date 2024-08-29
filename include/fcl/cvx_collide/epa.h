//
// Created by wei on 22-6-5.
//

#pragma once

#include "epa_polytope.h"
#include "epa_simplex2polytope.h"
#include "gjk.h"

namespace fcl {
namespace cvx_collide {

enum class EPA_Status { Failed, OK, Touching, IterationLimit, MallocFailed };

template <typename T>
class EPA {
 private:
  const std::size_t max_n_faces_;
  const std::size_t max_iterations_;
  const T tolerance_;

 public:
  explicit EPA(std::size_t max_n_faces = 127, std::size_t max_iterations = 255,
               T tolerance = 1e-6)
      : max_n_faces_(max_n_faces),
        max_iterations_(max_iterations),
        tolerance_(tolerance){};

  // Run EPA evaluation
  // The simplex (or gjk.simplex()) must contain the origin
  // The p_on_shape_0 and p_on_shape_1 are BOTH expressed IN SHAPE_0 frame
  // depth_if_penetration is POSITIVE
  EPA_Status Evaluate(const GJKSimplex<T>& simplex,
                      const MinkowskiDiff<T>& shape, T* depth_if_penetration,
                      Vector3<T>* p_on_shape_0, Vector3<T>* p_on_shape_1,
                      Polytope<T>* external_polytope_cache = nullptr) const;

  // Testing method
  EPA_Status Test_EvaluateFromInitializedPolytope(
      Polytope<T>& polytope, const MinkowskiDiff<T>& shape,
      T* depth_if_penetration, Vector3<T>* p_on_shape_0,
      Vector3<T>* p_on_shape_1) const;

  // Access of members
  std::size_t max_n_faces() const { return max_n_faces_; }
  std::size_t max_iterations() const { return max_iterations_; }
  T tolerance() const { return tolerance_; }

 private:
  // Internal evaluation interface
  EPA_Status evaluateFromInitializedPolytope(Polytope<T>& polytope,
                                             const MinkowskiDiff<T>& shape,
                                             T* depth_if_penetration,
                                             Vector3<T>* p_on_shape_0,
                                             Vector3<T>* p_on_shape_1) const;
  EPA_Status evaluateFromUnInitializedPolytope(const GJKSimplex<T>& simplex,
                                               Polytope<T>& polytope,
                                               const MinkowskiDiff<T>& shape,
                                               T* depth_if_penetration,
                                               Vector3<T>* p_on_shape_0,
                                               Vector3<T>* p_on_shape_1) const;

  // Determine the explore direction
  enum class FindNextDirectionStatus { OK, Failed, Converge };
  FindNextDirectionStatus findNextSupportDirection(
      const Polytope<T>& polytope, const MinkowskiDiff<T>& shape,
      PolytopeElementBase* nearest_feature, T tolerance, Vector3<T>& next_d,
      Vector3<T>& next_v, PolytopeFace<T>** start_face) const;

  // Get the min-distance point
  struct RawFeature {
    PolytopeElementType type;
    MinDistanceToSimplex<T> min_distance;
    MinkowskiDiffVertex<T> a_vertex, b_vertex, c_vertex;
  };
  static bool transformToRawFeature(const PolytopeElementBase* feature,
                                    RawFeature& raw_feature);
  bool checkTerminateCondition(const MinkowskiDiff<T>& shape,
                               const RawFeature& nearest_feature,
                               const Vector3<T>& d, const Vector3<T>& new_v,
                               T tolerance) const;
  void assignPenetrationPair(const MinkowskiDiff<T>& shape,
                             const MinkowskiDiffVertex<T>& candidate_next_v,
                             const RawFeature& nearest_feature,
                             T* depth_if_penetration, Vector3<T>* p_on_shape_0,
                             Vector3<T>* p_on_shape_1) const;
  bool assignPenetrationPairFromSegment(const MinkowskiDiff<T>& shape,
                                        const RawFeature& nearest_edge_feature,
                                        T* depth_if_penetration,
                                        Vector3<T>* p_on_shape_0,
                                        Vector3<T>* p_on_shape_1) const;
  bool assignPenetrationPairFromFace(const MinkowskiDiff<T>& shape,
                                     const RawFeature& nearest_face_feature,
                                     T* depth_if_penetration,
                                     Vector3<T>* p_on_shape_0,
                                     Vector3<T>* p_on_shape_1) const;
};

}  // namespace cvx_collide
}  // namespace fcl

#include "epa.hpp"
