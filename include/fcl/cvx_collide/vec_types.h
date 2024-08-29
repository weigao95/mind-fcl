//
// Created by Wei Gao on 2023/1/19.
//

#pragma once

#include <Eigen/Eigen>
#include <array>
#include <cassert>
#include <stack>

namespace fcl {
namespace cvx_collide {

template <typename T>
using Vector3 = Eigen::Matrix<T, 3, 1>;

template <typename T>
using Matrix3 = Eigen::Matrix<T, 3, 3>;

template <typename T>
using Transform3 = Eigen::Transform<T, 3, Eigen::Isometry>;

}  // namespace cvx_collide
}  // namespace fcl