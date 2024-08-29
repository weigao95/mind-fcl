#pragma once

namespace fcl {
namespace detail {

template <typename T>
GJKGeometryData<T> constructGJKGeometry(const ShapeBase<T>* shape) {
  // Check ptr type
  if (shape == nullptr) {
    GJKGeometryData<T> invalid;
    invalid.shape_type = GJKShapeType::Invalid;
    return invalid;
  }

  // Not none
  assert(shape != nullptr);
  switch (shape->getNodeType()) {
    case GEOM_TRIANGLE: {
      const auto* casted_ptr = static_cast<const TriangleP<T>*>(shape);
      return triangleToGJK(casted_ptr);
    }
    case GEOM_TETRAHEDRON: {
      const auto* casted_ptr = static_cast<const Tetrahedron<T>*>(shape);
      return tetrahedronToGJK(casted_ptr);
    }
    case GEOM_BOX: {
      const auto* casted_ptr = static_cast<const Box<T>*>(shape);
      return boxToGJK(*casted_ptr);
    }
    case GEOM_SPHERE: {
      const auto* casted_ptr = static_cast<const Sphere<T>*>(shape);
      return sphereToGJK(*casted_ptr);
    }
    case GEOM_ELLIPSOID: {
      const auto* casted_ptr = static_cast<const Ellipsoid<T>*>(shape);
      return ellipsoidToGJK(*casted_ptr);
    }
    case GEOM_CAPSULE: {
      const auto* casted_ptr = static_cast<const Capsule<T>*>(shape);
      return capsuleToGJK(*casted_ptr);
    }
    case GEOM_CONE: {
      const auto* casted_ptr = static_cast<const Cone<T>*>(shape);
      return coneToGJK(*casted_ptr);
    }
    case GEOM_CYLINDER: {
      const auto* casted_ptr = static_cast<const Cylinder<T>*>(shape);
      return cylinderToGJK(*casted_ptr);
    }
    case GEOM_CONVEX: {
      const auto* casted_ptr = static_cast<const Convex<T>*>(shape);
      return convexToGJK(casted_ptr);
    }
    default: {
      GJKGeometryData<T> gjk_geometry;
      gjk_geometry.shape_type = GJKShapeType::Invalid;
      return gjk_geometry;
    }
  }
}

template <typename T>
Vector3<T> computeSupportExceptSweptVolume(
    const GJKGeometryData<T>& gjk_geometry, const Vector3<T>& direction) {
  switch (gjk_geometry.shape_type) {
    case GJKShapeType::Triangle:
      return triangleSupport(gjk_geometry, direction);
    case GJKShapeType::Tetrahedron:
      return tetrahedronSupport(gjk_geometry, direction);
    case GJKShapeType::Convex:
      return convexSupport(gjk_geometry, direction);
    case GJKShapeType::Box:
      return boxSupport(gjk_geometry, direction);
    case GJKShapeType::Sphere:
      return sphereSupport(gjk_geometry, direction);
    case GJKShapeType::Ellipsoid:
      return ellipsoidSupport(gjk_geometry, direction);
    case GJKShapeType::Capsule:
      return capsuleSupport(gjk_geometry, direction);
    case GJKShapeType::Cone:
      return coneSupport(gjk_geometry, direction);
    case GJKShapeType::Cylinder:
      return cylinderSupport(gjk_geometry, direction);
    default:
      return Vector3<T>::Zero();
  }
}

template <typename T>
Vector3<T> computeSupport(const GJKGeometryData<T>& gjk_geometry,
                          const Vector3<T>& direction) {
  switch (gjk_geometry.shape_type) {
    case GJKShapeType::ShapeSweptVolume: {
      // Internal must be an swept volume
      const auto* swept_shape =
          static_cast<const GJKGeometryData<T>*>(gjk_geometry.user_ptr);
      auto point = computeSupportExceptSweptVolume(*swept_shape, direction);
      if (direction.dot(gjk_geometry.vec3_data) > 0) {
        return point + gjk_geometry.vec3_data;
      } else {
        return point;
      }
    }
    default: {
      return computeSupportExceptSweptVolume(gjk_geometry, direction);
    }
  }
}

template <typename T>
Vector3<T> computeInteriorExceptSweptVolume(
    const GJKGeometryData<T>& gjk_geometry) {
  switch (gjk_geometry.shape_type) {
    case GJKShapeType::Triangle: {
      return gjk_geometry.vec3_data;
    }
    case GJKShapeType::Tetrahedron: {
      return gjk_geometry.vec3_data;
    }
    case GJKShapeType::Convex: {
      const auto* convex = static_cast<const Convex<T>*>(gjk_geometry.user_ptr);
      return convex->getInteriorPoint();
    }
    default: {
      return basicGeometryInterior(gjk_geometry);
    }
  }
}

template <typename T>
Vector3<T> computeInterior(const GJKGeometryData<T>& gjk_geometry) {
  switch (gjk_geometry.shape_type) {
    case GJKShapeType::ShapeSweptVolume: {
      const auto* swept_shape =
          static_cast<const GJKGeometryData<T>*>(gjk_geometry.user_ptr);
      constexpr T half_displacement = 0.5;
      return computeInteriorExceptSweptVolume<T>(*swept_shape) +
             half_displacement * gjk_geometry.vec3_data;
    }
    default: {
      return computeInteriorExceptSweptVolume(gjk_geometry);
    }
  }
}

}  // namespace detail
}  // namespace fcl