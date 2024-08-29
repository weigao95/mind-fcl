//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/narrowphase/detail/ccd/ccd_typedef.h"
#include "fcl/narrowphase/detail/gjk_solver_cvx.h"

namespace fcl {
namespace detail {

template <typename S>
struct TranslationalCollisionGJK {
  // Binary interface
  static bool CheckSweptVolumeCollisionBinary(
      const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
      const Vector3<S>& shape1_displacement_in_tf1,
      const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
      int max_iteration = 128, S tolerance = 1e-6);
  static bool CheckSweptVolumeCollisionBinary(
      const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& shape1_displacement,
      const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
      int max_iteration = 128, S tolerance = 1e-6);

  // Return one sample of time-of-collision, which might NOT be the first
  static std::pair<bool, S> CheckSweptVolumeCollisionOneTocSample(
      const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
      const Vector3<S>& shape1_displacement_in_tf1,
      const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
      int max_iteration = 128, S tolerance = 1e-6);
  static std::pair<bool, S> CheckSweptVolumeCollisionOneTocSample(
      const GJKGeometryData<S>& shape1, const Transform3<S>& tf1,
      const TranslationalDisplacement<S>& shape1_displacement,
      const GJKGeometryData<S>& shape2, const Transform3<S>& tf2,
      int max_iteration = 128, S tolerance = 1e-6);

 private:
  static void setupMinkowskiDiff(const GJKGeometryData<S>& swept_volume1,
                                 const Transform3<S>& tf1,
                                 const GJKGeometryData<S>& shape2,
                                 const Transform3<S>& tf2,
                                 MinkowskiDiff<S>& shape);
  static S computeOneTocSampleForIntersection(
      const Vector3<S>& shape1_displacement_in_tf1,
      const typename detail::MPR<S>::IntersectData& intersect_data);
};

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/gjk_ccd-inl.h"