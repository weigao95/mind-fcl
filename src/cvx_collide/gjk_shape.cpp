#include "fcl/cvx_collide/gjk_shape.h"

namespace fcl {
namespace cvx_collide {

template struct GJKGeometryData<float>;
template struct GJKGeometryData<double>;

template Vector3<float> basicGeometrySupport(
    const GJKGeometryData<float>& gjk_geometry,
    const Vector3<float>& direction);
template Vector3<double> basicGeometrySupport(
    const GJKGeometryData<double>& gjk_geometry,
    const Vector3<double>& direction);
template Vector3<float> basicGeometryInterior(
    const GJKGeometryData<float>& gjk_geometry);
template Vector3<double> basicGeometryInterior(
    const GJKGeometryData<double>& gjk_geometry);

}  // namespace cvx_collide
}  // namespace fcl