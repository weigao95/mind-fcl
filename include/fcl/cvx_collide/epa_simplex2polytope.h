//
// Created by wei on 22-6-5.
//

#pragma once

#include "epa_polytope.h"
#include "gjk.h"

namespace fcl {
namespace cvx_collide {

enum class Simplex2PolytopeStatus { OK, Touching, Failed };

template <typename T>
Simplex2PolytopeStatus simplexToPolytope(const GJKSimplex<T>& gjk_simplex,
                                         const MinkowskiDiff<T>& shape,
                                         Polytope<T>& polytope,
                                         Vector3<T>* p0_if_touching,
                                         Vector3<T>* p1_if_touching,
                                         T touching_threshold = T(1e-10));

template <typename T>
Simplex2PolytopeStatus simplexToPolytope4(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiffVertex<T>& d,
    const MinkowskiDiff<T>& shape, Polytope<T>& polytope,
    Vector3<T>* p0_if_touching, Vector3<T>* p1_if_touching,
    T touching_threshold = T(1e-10));

template <typename T>
Simplex2PolytopeStatus simplexToPolytope3(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiff<T>& shape,
    Polytope<T>& polytope, Vector3<T>* p0_if_touching,
    Vector3<T>* p1_if_touching, T touching_threshold = T(1e-10));

template <typename T>
Simplex2PolytopeStatus simplexToPolytope2(const MinkowskiDiffVertex<T>& a,
                                          const MinkowskiDiffVertex<T>& b,
                                          const MinkowskiDiff<T>& shape,
                                          Polytope<T>& polytope,
                                          Vector3<T>* p0_if_touching,
                                          Vector3<T>* p1_if_touching,
                                          T touching_threshold = T(1e-10));

template <typename T>
Simplex2PolytopeStatus formNewTetrahedronPolytope(
    const MinkowskiDiffVertex<T>& a, const MinkowskiDiffVertex<T>& b,
    const MinkowskiDiffVertex<T>& c, const MinkowskiDiffVertex<T>& d,
    Polytope<T>& polytope);

}  // namespace cvx_collide
}  // namespace fcl

#include "epa_simplex2polytope.hpp"
