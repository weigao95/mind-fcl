#pragma once

namespace fcl {
namespace cvx_collide {

//==============================================================================
template <typename T>
MinkowskiDiff<T>::MinkowskiDiff()
    : support_function(basicGeometrySupport<T>),
      interior_function(basicGeometryInterior<T>) {}

//==============================================================================
template <typename T>
Vector3<T> MinkowskiDiff<T>::support0(const Vector3<T>& d) const {
  return support_function(shapes[0], d);
}

//==============================================================================
template <typename T>
Vector3<T> MinkowskiDiff<T>::support1(const Vector3<T>& d) const {
#ifdef WIN32
  // Convert direction to 1
  const T direction_in_1_x =
      toshape1(0, 0) * d[0] + toshape1(0, 1) * d[1] + toshape1(0, 2) * d[2];
  const T direction_in_1_y =
      toshape1(1, 0) * d[0] + toshape1(1, 1) * d[1] + toshape1(1, 2) * d[2];
  const T direction_in_1_z =
      toshape1(2, 0) * d[0] + toshape1(2, 1) * d[1] + toshape1(2, 2) * d[2];
  const Vector3<T> direction_in_1(direction_in_1_x, direction_in_1_y,
                                  direction_in_1_z);
  const Vector3<T> support_in_1 = support_function(shapes[1], direction_in_1);

  // Convert back to 0
  const auto& to_shape0_matrix = toshape0.matrix();
  const T support_in_0_x = to_shape0_matrix(0, 0) * support_in_1[0] +
                           to_shape0_matrix(0, 1) * support_in_1[1] +
                           to_shape0_matrix(0, 2) * support_in_1[2] +
                           to_shape0_matrix(0, 3);
  const T support_in_0_y = to_shape0_matrix(1, 0) * support_in_1[0] +
                           to_shape0_matrix(1, 1) * support_in_1[1] +
                           to_shape0_matrix(1, 2) * support_in_1[2] +
                           to_shape0_matrix(1, 3);
  const T support_in_0_z = to_shape0_matrix(2, 0) * support_in_1[0] +
                           to_shape0_matrix(2, 1) * support_in_1[1] +
                           to_shape0_matrix(2, 2) * support_in_1[2] +
                           to_shape0_matrix(2, 3);
  return Vector3<T>(support_in_0_x, support_in_0_y, support_in_0_z);
#else
  const Vector3<T> dir_in_1 = toshape1 * d;
  return toshape0 * support_function(shapes[1], dir_in_1);
#endif
}

//==============================================================================
template <typename T>
Vector3<T> MinkowskiDiff<T>::support(const Vector3<T>& d) const {
  return support0(d) - support1(-d);
}

//==============================================================================
template <typename T>
MinkowskiDiffVertex<T> MinkowskiDiff<T>::supportVertex(
    const Vector3<T>& d) const {
  MinkowskiDiffVertex<T> v;
  v.direction = d;
  v.vertex = support0(d) - support1(-d);
  return v;
}

//==============================================================================
template <typename T>
Vector3<T> MinkowskiDiff<T>::support(const Vector3<T>& d, size_t index) const {
  if (index)
    return support1(d);
  else
    return support0(d);
}

//==============================================================================
template <typename T>
Vector3<T> MinkowskiDiff<T>::interior() const {
  Vector3<T> shape_0_interior = interior_function(shapes[0]);
  Vector3<T> shape_1_interior = interior_function(shapes[1]);
  return shape_0_interior - toshape0 * shape_1_interior;
}

}  // namespace cvx_collide
}  // namespace fcl