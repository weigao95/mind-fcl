#pragma once

namespace fcl {
namespace cvx_collide {

// Interface for evaluation
template <typename T>
Vector3<T> basicGeometrySupport(const GJKGeometryData<T>& gjk_geometry,
                                const Vector3<T>& direction) {
  switch (gjk_geometry.shape_type) {
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
Vector3<T> basicGeometryInterior(const GJKGeometryData<T>& gjk_geometry) {
  (void)(gjk_geometry);
  return Vector3<T>::Zero();
}

// For box type
template <typename T>
GJKGeometryData<T> constructGJKBox(const Vector3<T>& box_size_xyz) {
  return GJKGeometryData<T>(GJKShapeType::Box, box_size_xyz);
}

template <typename T>
Vector3<T> boxSupport(const GJKGeometryData<T>& geometry_data,
                      const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Box);
  const Vector3<T>& side = geometry_data.vec3_data;
  return Vector3<T>((dir[0] > 0) ? (side[0] / 2) : (-side[0] / 2),
                    (dir[1] > 0) ? (side[1] / 2) : (-side[1] / 2),
                    (dir[2] > 0) ? (side[2] / 2) : (-side[2] / 2));
}

// For sphere type
template <typename T>
GJKGeometryData<T> constructGJKSphere(T sphere_radius) {
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Sphere;
  Vector3<T>& vec3_data = geometry_data.vec3_data;
  vec3_data[0] = sphere_radius;
  return geometry_data;
}

template <typename T>
Vector3<T> sphereSupport(const GJKGeometryData<T>& geometry_data,
                         const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Sphere);
  const Vector3<T>& vec3_data = geometry_data.vec3_data;
  const T radius = vec3_data[0];
  return dir * radius;
}

// For ellipsoid type
template <typename T>
GJKGeometryData<T> constructGJKEllipsoid(const Vector3<T>& radii) {
  return GJKGeometryData<T>(GJKShapeType::Ellipsoid, radii);
}

template <typename T>
Vector3<T> ellipsoidSupport(const GJKGeometryData<T>& geometry_data,
                            const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Ellipsoid);
  const Vector3<T>& radii = geometry_data.vec3_data;
  const T a2 = radii[0] * radii[0];
  const T b2 = radii[1] * radii[1];
  const T c2 = radii[2] * radii[2];
  const Vector3<T> v(a2 * dir[0], b2 * dir[1], c2 * dir[2]);
  const T d = std::sqrt(v.dot(dir));
  return v / d;
}

// For capsule type
template <typename T>
GJKGeometryData<T> constructGJKCapsule(T radius, T length_z) {
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Capsule;
  Vector3<T>& vec3_data = geometry_data.vec3_data;
  vec3_data[0] = radius;
  vec3_data[1] = length_z;
  return geometry_data;
}

template <typename T>
Vector3<T> capsuleSupport(const GJKGeometryData<T>& geometry_data,
                          const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Capsule);
  const Vector3<T>& vec3_data = geometry_data.vec3_data;
  const T radius = vec3_data[0];
  const T lz = vec3_data[1];
  const T half_h = lz * 0.5;
  Vector3<T> pos1(0, 0, half_h);
  Vector3<T> pos2(0, 0, -half_h);
  Vector3<T> v = dir * radius;
  pos1 += v;
  pos2 += v;
  if (dir.dot(pos1) > dir.dot(pos2))
    return pos1;
  else
    return pos2;
}

// For cone type
template <typename T>
GJKGeometryData<T> constructGJKCone(T radius, T length_z) {
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Cone;
  Vector3<T>& vec3_data = geometry_data.vec3_data;
  vec3_data[0] = radius;
  vec3_data[1] = length_z;
  return geometry_data;
}

template <typename T>
Vector3<T> coneSupport(const GJKGeometryData<T>& geometry_data,
                       const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Cone);
  const Vector3<T>& vec3_data = geometry_data.vec3_data;
  const T radius = vec3_data[0];
  const T lz = vec3_data[1];

  T zdist = dir[0] * dir[0] + dir[1] * dir[1];
  T len = zdist + dir[2] * dir[2];
  zdist = std::sqrt(zdist);
  len = std::sqrt(len);
  const T half_h = lz * 0.5;
  const T sin_a = radius / std::sqrt(radius * radius + 4 * half_h * half_h);

  if (dir[2] > len * sin_a) {
    return Vector3<T>(0, 0, half_h);
  } else if (zdist > 0) {
    const T rad = radius / zdist;
    return Vector3<T>(rad * dir[0], rad * dir[1], -half_h);
  } else {
    return Vector3<T>(0, 0, -half_h);
  }
}

// For cylinder type
template <typename T>
GJKGeometryData<T> constructGJKCylinder(T radius, T length_z) {
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Cylinder;
  Vector3<T>& vec3_data = geometry_data.vec3_data;
  vec3_data[0] = radius;
  vec3_data[1] = length_z;
  return geometry_data;
}

template <typename T>
Vector3<T> cylinderSupport(const GJKGeometryData<T>& geometry_data,
                           const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Cylinder);
  const Vector3<T>& vec3_data = geometry_data.vec3_data;
  const T radius = vec3_data[0];
  const T lz = vec3_data[1];
  const T zdist = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
  const T half_h = lz * 0.5;
  if (zdist == 0.0) {
    return Vector3<T>(0, 0, (dir[2] > 0) ? half_h : -half_h);
  } else {
    const T d = radius / zdist;
    return Vector3<T>(d * dir[0], d * dir[1], (dir[2] > 0) ? half_h : -half_h);
  }
}

}  // namespace cvx_collide
}  // namespace fcl