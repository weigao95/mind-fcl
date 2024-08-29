//
// Created by wei on 2/1/22.
//

#pragma once

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/cylinder.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/geometry/shape/convex.h"

namespace fcl {
namespace test {

template <typename BV>
::fcl::BVHModel<BV>* makeCubeBVH(typename BV::S size_x, typename BV::S size_y,
                                 typename BV::S size_z);

template <typename T>
void meshEllipsoidTriangles(T size_x, T size_y, T size_z,
                            std::vector<std::array<T, 3>>* points,
                            std::vector<std::array<int, 3>>* triangles);
template <typename T>
void meshEllipsoidUV(T size_x, T size_y, T size_z, int n_slices, int n_stacks,
                     std::vector<std::array<T, 3>>* points,
                     std::vector<std::array<int, 3>>* triangles);

template <typename BV>
::fcl::BVHModel<BV>* makeEllipsoidBVH(typename BV::S size_x,
                                      typename BV::S size_y,
                                      typename BV::S size_z);
template <typename T>
::fcl::Convex<T> makeEllipsoidConvex(T size_x, T size_y, T size_z);
template <typename T>
::fcl::Convex<T> makeEllipsoidConvexUV(T size_x, T size_y, T size_z,
                                       int n_slices, int n_stacks);

template <typename BV>
::fcl::BVHModel<BV>* makeSphereBVH(typename BV::S radius) {
  return makeEllipsoidBVH<BV>(2.0 * radius, 2.0 * radius, 2.0 * radius);
}

template <typename T>
::fcl::Convex<T> makeSphereConvex(T radius) {
  return makeEllipsoidConvex<T>(2.0 * radius, 2.0 * radius, 2.0 * radius);
}
template <typename T>
::fcl::Convex<T> makeSphereConvexUV(T radius, int n_slices, int n_stacks) {
  return makeEllipsoidConvexUV<T>(2.0 * radius, 2.0 * radius, 2.0 * radius,
                                  n_slices, n_stacks);
}

}  // namespace test
}  // namespace fcl

// template implementation in header
#include "create_primitive_mesh-inl.h"
