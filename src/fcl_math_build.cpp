//
// Created by mech-mind_gw on 4/23/2023.
//

#include "fcl/math/bv/AABB.h"
#include "fcl/math/bv/OBB.h"
#include "fcl/math/bv/OBBRSS.h"
#include "fcl/math/bv/RSS.h"
#include "fcl/math/bv/kDOP.h"
#include "fcl/math/bv/kIOS.h"
#include "fcl/math/bv/utility.h"
#include "fcl/math/constants.h"
#include "fcl/math/geometry.h"
#include "fcl/math/math_simd_details.h"
#include "fcl/math/rng.h"
#include "fcl/math/mesh_simplex.h"
#include "fcl/math/variance3.h"

namespace fcl {

// math/bv/*.cpp
//==============================================================================
template class AABB<float>;

// skip kDOP/kIOS

template class OBB<float>;

template void computeVertices(const OBB<float> &b, Vector3<float> vertices[8]);

template OBB<float> merge_largedist(const OBB<float> &b1, const OBB<float> &b2);

template OBB<float> merge_smalldist(const OBB<float> &b1, const OBB<float> &b2);

template class OBBRSS<float>;

template OBBRSS<float> translate(const OBBRSS<float> &bv,
                                 const Vector3<float> &t);

template class RSS<float>;

template void clipToRange(float &val, float a, float b);

template void segCoords(float &t, float &u, float a, float b, float A_dot_B,
                        float A_dot_T, float B_dot_T);

template bool inVoronoi(float a, float b, float Anorm_dot_B, float Anorm_dot_T,
                        float A_dot_B, float A_dot_T, float B_dot_T);

template float rectDistance(const Matrix3<float> &Rab,
                            const Vector3<float> &Tab, const float a[2],
                            const float b[2], Vector3<float> *P,
                            Vector3<float> *Q);

template float rectDistance(const Transform3<float> &tfab, const float a[2],
                            const float b[2], Vector3<float> *P,
                            Vector3<float> *Q);

template RSS<float> translate(const RSS<float> &bv, const Vector3<float> &t);

// math/bv/*.cpp
//==============================================================================
template class AABB<double>;

// skip kDOP/kIOS

template class OBB<double>;

template void computeVertices(const OBB<double> &b,
                              Vector3<double> vertices[8]);

template OBB<double> merge_largedist(const OBB<double> &b1,
                                     const OBB<double> &b2);

template OBB<double> merge_smalldist(const OBB<double> &b1,
                                     const OBB<double> &b2);

template class OBBRSS<double>;

template OBBRSS<double> translate(const OBBRSS<double> &bv,
                                  const Vector3<double> &t);

template class RSS<double>;

template void clipToRange(double &val, double a, double b);

template void segCoords(double &t, double &u, double a, double b,
                        double A_dot_B, double A_dot_T, double B_dot_T);

template bool inVoronoi(double a, double b, double Anorm_dot_B,
                        double Anorm_dot_T, double A_dot_B, double A_dot_T,
                        double B_dot_T);

template double rectDistance(const Matrix3<double> &Rab,
                             const Vector3<double> &Tab, const double a[2],
                             const double b[2], Vector3<double> *P,
                             Vector3<double> *Q);

template double rectDistance(const Transform3<double> &tfab, const double a[2],
                             const double b[2], Vector3<double> *P,
                             Vector3<double> *Q);

template RSS<double> translate(const RSS<double> &bv, const Vector3<double> &t);

} // namespace fcl