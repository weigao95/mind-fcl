#pragma once

namespace fcl {
namespace detail {

template <typename Shape>
void initializeShapeFixedOrientationBoxTranslationalCCD(
    const Shape& shape1, const Transform3<typename Shape::S>& tf_shape1,
    const TranslationalDisplacement<typename Shape::S>& displacement_shape1,
    const Transform3<typename Shape::S>& tf_another_2,
    // output
    FixedOrientationBoxPairTranslationalCCD<typename Shape::S>& disjoint,
    AABB<typename Shape::S>& shape_local_AABB) {
  using S = typename Shape::S;
  // Compute the bv for shape
  OBB<S> shape_obb_world;
  computeBV(shape1, tf_shape1, shape_obb_world);

  // Convert OBB to local AABB and a tf_AABB on that AABB
  shape_local_AABB.max_ = shape_obb_world.extent;
  shape_local_AABB.min_ = -shape_obb_world.extent;

  // Make tf_AABB frame
  Transform3<S> tf_shape_AABB;
  tf_shape_AABB.setIdentity();
  tf_shape_AABB.linear().matrix() = shape_obb_world.axis;
  tf_shape_AABB.translation() = shape_obb_world.To;

  TranslationalDisplacement<S> displacement_box;
  displacement_box.scalar_displacement =
      displacement_shape1.scalar_displacement;
  displacement_box.unit_axis_in_shape1 =
      tf_shape_AABB.linear().transpose() *
      (tf_shape1.linear() * displacement_shape1.unit_axis_in_shape1);
  disjoint.Initialize(tf_shape_AABB, tf_another_2, displacement_box);
}

}  // namespace detail
}  // namespace fcl