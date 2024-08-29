#pragma once
#include "gjk_shape.h"

namespace fcl {
namespace cvx_collide {

/// @brief The vertex/direction to obtain that vertex on Mink
template <typename T>
struct MinkowskiDiffVertex {
  Vector3<T> vertex;
  Vector3<T> direction;
};

/// @brief Minkowski difference class of two shapes
template <typename T>
struct MinkowskiDiff {
  /// @brief points to two shapes
  GJKGeometryData<T> shapes[2];

  /// Function for support and interior
  /// Default value are provided, can be over-written if necessary
  SupportFunction<T> support_function;
  InteriorFunction<T> interior_function;

  /// @brief rotation from shape0 to shape1
  Matrix3<T> toshape1;

  /// @brief transform from shape1 to shape0
  Transform3<T> toshape0;

  MinkowskiDiff();

  /// @brief support function for shape0
  Vector3<T> support0(const Vector3<T>& d) const;

  /// @brief support function for shape1
  Vector3<T> support1(const Vector3<T>& d) const;

  /// @brief support function for the pair of shapes
  Vector3<T> support(const Vector3<T>& d) const;
  MinkowskiDiffVertex<T> supportVertex(const Vector3<T>& d) const;

  /// @brief support function for the d-th shape (d = 0 or 1)
  Vector3<T> support(const Vector3<T>& d, size_t index) const;

  /// @brief Compute an interior point inside the MinkowskiDiff
  Vector3<T> interior() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace cvx_collide
}  // namespace fcl

#include "minkowski_diff.hpp"
