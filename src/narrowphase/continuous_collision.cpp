//
// Created by wei on 24-6-18.
//
#include "fcl/narrowphase/continuous_collision.h"

namespace fcl {

template struct ContinuousCollisionContact<float>;
template struct ContinuousCollisionContact<double>;

template struct ContinuousCollisionRequest<float>;
template struct ContinuousCollisionRequest<double>;

template struct ContinuousCollisionResult<float>;
template struct ContinuousCollisionResult<double>;

template void translational_ccd(
    const CollisionGeometry<float>* o1, const Transform3<float>& tf1,
    const TranslationalDisplacement<float>& o1_displacement,
    const CollisionGeometry<float>* o2, const Transform3<float>& tf2,
    const ContinuousCollisionRequest<float>& request,
    ContinuousCollisionResult<float>& result);
template void translational_ccd(
    const CollisionGeometry<double>* o1, const Transform3<double>& tf1,
    const TranslationalDisplacement<double>& o1_displacement,
    const CollisionGeometry<double>* o2, const Transform3<double>& tf2,
    const ContinuousCollisionRequest<double>& request,
    ContinuousCollisionResult<double>& result);

} // namespace fcl