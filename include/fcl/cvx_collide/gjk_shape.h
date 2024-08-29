//
// Created by wei on 22-11-11.
//

#pragma once
#include "vec_types.h"

namespace fcl {
namespace cvx_collide {

enum class GJKShapeType {
  // Basic geometry types
  Box = 0,
  Sphere = 1,
  Ellipsoid = 2,
  Capsule = 3,
  Cone = 4,
  Cylinder = 5,

  // Compound types via ptr access
  Convex = 100,
  Triangle = 101,
  Tetrahedron = 102,
  ShapeSweptVolume = 103,

  // Error
  Invalid = 1000,
};

// The data type to hold all the geometry
// The liveness of this class is very short, usually
// not more than the GJK::Evaluate class, or no more
// than the user_ptr which is the compound shapes
template <typename T>
struct GJKGeometryData {
  GJKShapeType shape_type{GJKShapeType::Invalid};
  Vector3<T> vec3_data;
  const void* user_ptr{nullptr};

  // Construct from vector
  explicit GJKGeometryData() = default;
  GJKGeometryData(GJKShapeType shape_type, Vector3<T> vec_data)
      : shape_type(shape_type), vec3_data(std::move(vec_data)) {}
};

// Function for support and interior
template <typename T>
using SupportFunction =
    std::function<Vector3<T>(const GJKGeometryData<T>&, const Vector3<T>&)>;
template <typename T>
using InteriorFunction = std::function<Vector3<T>(const GJKGeometryData<T>&)>;

// Evaluate the support and interior
template <typename T>
Vector3<T> basicGeometrySupport(const GJKGeometryData<T>& gjk_geometry,
                                const Vector3<T>& direction);
template <typename T>
Vector3<T> basicGeometryInterior(const GJKGeometryData<T>& gjk_geometry);

// For different types
template <typename T>
GJKGeometryData<T> constructGJKBox(const Vector3<T>& box_size_xyz);
template <typename T>
Vector3<T> boxSupport(const GJKGeometryData<T>& geometry_data,
                      const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> constructGJKSphere(T sphere_radius);
template <typename T>
Vector3<T> sphereSupport(const GJKGeometryData<T>& geometry_data,
                         const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> constructGJKEllipsoid(const Vector3<T>& radii);
template <typename T>
Vector3<T> ellipsoidSupport(const GJKGeometryData<T>& geometry_data,
                            const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> constructGJKCapsule(T radius, T length_z);
template <typename T>
Vector3<T> capsuleSupport(const GJKGeometryData<T>& geometry_data,
                          const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> constructGJKCone(T radius, T length_z);
template <typename T>
Vector3<T> coneSupport(const GJKGeometryData<T>& geometry_data,
                       const Vector3<T>& direction);

template <typename T>
GJKGeometryData<T> constructGJKCylinder(T radius, T length_z);
template <typename T>
Vector3<T> cylinderSupport(const GJKGeometryData<T>& geometry_data,
                           const Vector3<T>& direction);

}  // namespace cvx_collide
}  // namespace fcl

#include "gjk_shape.hpp"
