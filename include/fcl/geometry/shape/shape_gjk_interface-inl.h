#pragma once

namespace fcl {
namespace detail {

// New shape types
template <typename T>
GJKGeometryData<T> triangleToGJK(const TriangleP<T>* triangle) {
  assert(triangle != nullptr);
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Triangle;
  geometry_data.user_ptr = triangle;
  geometry_data.vec3_data = (triangle->a + triangle->b + triangle->c) / T(3.0);
  return geometry_data;
}

template <typename T>
Vector3<T> triangleSupport(const GJKGeometryData<T>& geometry_data,
                           const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Triangle);
  const auto* triangle =
      static_cast<const TriangleP<T>*>(geometry_data.user_ptr);
  const Vector3<T>& a = triangle->a;
  const Vector3<T>& b = triangle->b;
  const Vector3<T>& c = triangle->c;

  const T dota = dir.dot(a);
  const T dotb = dir.dot(b);
  const T dotc = dir.dot(c);
  if (dota > dotb) {
    if (dotc > dota)
      return c;
    else
      return a;
  } else {
    if (dotc > dotb)
      return c;
    else
      return b;
  }
}

template <typename T>
GJKGeometryData<T> tetrahedronToGJK(const Tetrahedron<T>* tetrahedron) {
  assert(tetrahedron != nullptr);
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Tetrahedron;
  geometry_data.user_ptr = tetrahedron;
  const auto& vertices = tetrahedron->vertices;
  geometry_data.vec3_data =
      (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / T(4.0);
  return geometry_data;
}

template <typename T>
Vector3<T> tetrahedronSupport(const GJKGeometryData<T>& geometry_data,
                              const Vector3<T>& dir) {
  assert(geometry_data.shape_type == GJKShapeType::Tetrahedron);
  const auto* tetrahedron =
      static_cast<const Tetrahedron<T>*>(geometry_data.user_ptr);
  const auto& vertices = tetrahedron->vertices;

  T max_dot = -std::numeric_limits<T>::infinity();
  int max_dot_index = -1;
  for (int i = 0; i < 4; i++) {
    const T dot_i = dir.dot(vertices[i]);
    if (dot_i > max_dot) {
      max_dot_index = i;
      max_dot = dot_i;
    }
  }

  return vertices[max_dot_index];
}

// For convex type
template <typename T>
GJKGeometryData<T> convexToGJK(const Convex<T>* convex) {
  assert(convex != nullptr);
  GJKGeometryData<T> geometry_data;
  geometry_data.shape_type = GJKShapeType::Convex;
  geometry_data.user_ptr = convex;
  return geometry_data;
}

template <typename T>
Vector3<T> convexSupport(const GJKGeometryData<T>& geometry_data,
                         const Vector3<T>& direction) {
  assert(geometry_data.shape_type == GJKShapeType::Convex);
  const auto* convex = static_cast<const Convex<T>*>(geometry_data.user_ptr);
  return convex->findExtremeVertex(direction);
}

// Construct existing shapes
template <typename T>
GJKGeometryData<T> boxToGJK(const Box<T>& box) {
  return cvx_collide::constructGJKBox(box.side);
}

template <typename T>
GJKGeometryData<T> sphereToGJK(const Sphere<T>& sphere) {
  return cvx_collide::constructGJKSphere(sphere.radius);
}

template <typename T>
GJKGeometryData<T> ellipsoidToGJK(const Ellipsoid<T>& ellipsoid) {
  return cvx_collide::constructGJKEllipsoid(ellipsoid.radii);
}

template <typename T>
GJKGeometryData<T> capsuleToGJK(const Capsule<T>& capsule) {
  return cvx_collide::constructGJKCapsule(capsule.radius, capsule.lz);
}

template <typename T>
GJKGeometryData<T> coneToGJK(const Cone<T>& cone) {
  return cvx_collide::constructGJKCone(cone.radius, cone.lz);
}

template <typename T>
GJKGeometryData<T> cylinderToGJK(const Cylinder<T>& cylinder) {
  return cvx_collide::constructGJKCylinder(cylinder.radius, cylinder.lz);
}

/// templated shape2gjk
template <typename S, typename Shape>
struct ShapeToGJK_Impl {
  static GJKGeometryData<S> run(const Shape&) {
    auto invalid = GJKGeometryData<S>();
    invalid.shape_type = GJKShapeType::Invalid;
    return invalid;
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Box<S>> {
  static GJKGeometryData<S> run(const Box<S>& shape) {
    return boxToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Sphere<S>> {
  static GJKGeometryData<S> run(const Sphere<S>& shape) {
    return sphereToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Ellipsoid<S>> {
  static GJKGeometryData<S> run(const Ellipsoid<S>& shape) {
    return ellipsoidToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Capsule<S>> {
  static GJKGeometryData<S> run(const Capsule<S>& shape) {
    return capsuleToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Cone<S>> {
  static GJKGeometryData<S> run(const Cone<S>& shape) {
    return coneToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Cylinder<S>> {
  static GJKGeometryData<S> run(const Cylinder<S>& shape) {
    return cylinderToGJK<S>(shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, Convex<S>> {
  static GJKGeometryData<S> run(const Convex<S>& shape) {
    return convexToGJK<S>(&shape);
  }
};

//==============================================================================
template <typename S>
struct ShapeToGJK_Impl<S, TriangleP<S>> {
  static GJKGeometryData<S> run(const TriangleP<S>& shape) {
    return triangleToGJK<S>(shape);
  }
};

//==============================================================================
template <typename Shape>
GJKGeometryData<typename Shape::S> shapeToGJK(const Shape& shape) {
  return ShapeToGJK_Impl<typename Shape::S, Shape>::run(shape);
}

}  // namespace detail
}  // namespace fcl