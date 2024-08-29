//
// Created by mech-mind_gw on 3/31/2023.
//

#pragma once

#include "fcl/narrowphase/detail/gjk_solver_cvx.h"

namespace fcl {
namespace detail {

template <typename S>
struct PenetrationObjectCallbackCache {
  // For shape
  GJKGeometryData<S> shape_gjk_data;

  // For bvh
  const MeshSimplex* bvh_simplexes{nullptr};
  const Vector3<S>* bvh_vertices{nullptr};
  TriangleP<S> bvh_triangle_instance;
  Tetrahedron<S> bvh_tetrahedron_instance;

  explicit PenetrationObjectCallbackCache()
      : bvh_triangle_instance(Vector3<S>::Zero(), Vector3<S>::Zero(),
                              Vector3<S>::Zero()),
        bvh_tetrahedron_instance(Vector3<S>::Zero(), Vector3<S>::Zero(),
                                 Vector3<S>::Zero(), Vector3<S>::Zero()) {
    shape_gjk_data.shape_type = GJKShapeType::Invalid;
  }
};

//==============================================================================
template <typename S>
void penetrationDistanceGetContactGJK(
    const Contact<S>& contact, const CollisionGeometry<S>* geom,
    const Transform3<S>& tf_geom, bool geom_is_1, GJKGeometryData<S>& gjk_geom,
    Transform3<S>& tf_gjk_geom, PenetrationObjectCallbackCache<S>& cache) {
  const bool is_o1 =
      (contact.o1 == contact.o2) ? geom_is_1 : (contact.o1 == geom);
  switch (geom->getObjectType()) {
    case OBJECT_TYPE::OT_GEOM: {
      if (cache.shape_gjk_data.shape_type == GJKShapeType::Invalid) {
        cache.shape_gjk_data =
            constructGJKGeometry(static_cast<const ShapeBase<S>*>(geom));
      }

      // Assign the result
      gjk_geom = cache.shape_gjk_data;
      tf_gjk_geom = tf_geom;
      assert(gjk_geom.shape_type != GJKShapeType::Invalid);
      return;
    }
    case OBJECT_TYPE::OT_OCTREE2:
    case OBJECT_TYPE::OT_HEIGHTMAP: {
      const AABB<S>& aabb_box = is_o1 ? contact.o1_bv : contact.o2_bv;
      const Vector3<S> box_side = aabb_box.max_ - aabb_box.min_;
      gjk_geom = ::fcl::cvx_collide::constructGJKBox(box_side);
      tf_gjk_geom = tf_geom;
      tf_gjk_geom.translation() += tf_geom.linear() * aabb_box.center();
      return;
    }
    case OBJECT_TYPE::OT_BVH: {
      // Init the cache if not yet
      if (!cache.bvh_simplexes) {
        const auto* bvh = static_cast<const fcl::BVHModelBase<S>*>(geom);
        cache.bvh_simplexes = bvh->getSimplex();
        cache.bvh_vertices = bvh->getVertices();
      }

      // Make the triangle
      const MeshSimplex& mesh_simplex =
          cache.bvh_simplexes[is_o1 ? contact.b1 : contact.b2];
      if (mesh_simplex.is_triangle()) {
        cache.bvh_triangle_instance.a = cache.bvh_vertices[mesh_simplex[0]];
        cache.bvh_triangle_instance.b = cache.bvh_vertices[mesh_simplex[1]];
        cache.bvh_triangle_instance.c = cache.bvh_vertices[mesh_simplex[2]];
        gjk_geom = constructGJKGeometry(&cache.bvh_triangle_instance);
      } else {
        for (int i = 0; i < 4; i++) {
          cache.bvh_tetrahedron_instance.vertices[i] =
              cache.bvh_vertices[mesh_simplex[i]];
        }
        gjk_geom = constructGJKGeometry(&cache.bvh_tetrahedron_instance);
      }

      // Pose of simplex is the same as tf_geom
      tf_gjk_geom = tf_geom;
      return;
    }
    default:
      assert(false);
  }
}

//==============================================================================
template <typename S>
void penetrationDepthSwapContactObjectPair(Contact<S>& contact) {
  using std::swap;
  swap(contact.o1, contact.o2);
  swap(contact.b1, contact.b2);
  swap(contact.o1_bv, contact.o2_bv);
  contact.normal *= -1;
}

//==============================================================================
template <typename S>
void computePenetrationMPR(
    const GJKGeometryData<S>& geom1, const Transform3<S>& tf_1,
    const GJKGeometryData<S>& geom2, const Transform3<S>& tf_2,
    const Vector3<S>& geom2_escape_direction_world_init, const MPR<S>& mpr,
    bool is_incremental, Vector3<S>& contact_position_out,
    Vector3<S>& contact_normal_out, S& penetration_depth_out) {
  if (geom1.shape_type == GJKShapeType::Invalid ||
      geom2.shape_type == GJKShapeType::Invalid) {
    contact_position_out.setZero();
    penetration_depth_out = S(-1);
    return;
  }

  // Construct MinkowskiDiff
  MinkowskiDiff<S> shape;
  shape.shapes[0] = geom1;
  shape.shapes[1] = geom2;
  shape.support_function = computeSupport<S>;
  shape.interior_function = computeInterior<S>;
  shape.toshape1.noalias() = tf_2.linear().transpose() * tf_1.linear();
  shape.toshape0 = tf_1.inverse(Eigen::Isometry) * tf_2;

  // Map direction to shape1 frame, in which the MinkowskiDiff lies
  Vector3<S> geom2_escape_direction =
      tf_1.linear().transpose() * geom2_escape_direction_world_init;

  // Invoke penetration
  if (is_incremental) {
    using IncrementalMinimumPenetrationData =
        typename ::fcl::cvx_collide::MPR<S>::IncrementalMinimumPenetrationData;
    using IncrementalPenetrationStatus =
        typename ::fcl::cvx_collide::MPR<S>::IncrementalPenetrationStatus;
    IncrementalMinimumPenetrationData penetration_data;
    IncrementalPenetrationStatus status =
        mpr.IncrementalMinimumPenetrationDistance(shape, geom2_escape_direction,
                                                  penetration_data);
    if (status == IncrementalPenetrationStatus::OK ||
        // At iteration limit, it is still a valid penetration
        // Although it might be not optimal
        status == IncrementalPenetrationStatus::IterationLimit) {
      penetration_depth_out = penetration_data.minimum_penetration;
      Vector3<S> contact_center =
          S(0.5) * (penetration_data.p0_in_shape0_frame +
                    penetration_data.p1_in_shape0_frame);
      contact_position_out = tf_1 * contact_center;
      contact_normal_out =
          tf_1.linear() * penetration_data.penetration_direction;
      return;
    } else {
      contact_position_out.setZero();
      penetration_depth_out = S(-1);
      contact_normal_out = geom2_escape_direction_world_init;
      return;
    }
  } else {
    // Directed penetration
    using DirectedPenetrationData =
        typename ::fcl::cvx_collide::MPR<S>::DirectedPenetrationData;
    using DirectedPenetrationStatus =
        typename ::fcl::cvx_collide::MPR<S>::DirectedPenetrationStatus;
    DirectedPenetrationData penetration_data;
    DirectedPenetrationStatus status = mpr.DirectedPenetration(
        shape, geom2_escape_direction, penetration_data);
    if (status == DirectedPenetrationStatus::OK) {
      penetration_depth_out = penetration_data.distance_on_direction;
      Vector3<S> contact_center =
          S(0.5) * (penetration_data.p0_in_shape0_frame +
                    penetration_data.p1_in_shape0_frame);
      contact_position_out = tf_1 * contact_center;
      contact_normal_out = geom2_escape_direction_world_init;
      return;
    } else {
      contact_position_out.setZero();
      penetration_depth_out = S(-1);
      contact_normal_out = geom2_escape_direction_world_init;
      return;
    }
  }
}

//==============================================================================
template <typename S>
std::size_t collisionPenetrationMPR(const CollisionGeometry<S>* o1,
                                    const Transform3<S>& tf1,
                                    const CollisionGeometry<S>* o2,
                                    const Transform3<S>& tf2,
                                    const GJKSolver<S>& solver,
                                    const CollisionRequest<S>& request,
                                    CollisionResult<S>& result) {
  // Direction request
  const auto& penetration_mode = request.penetrationMode();
  const Vector3<S>& shape2_escape_direction_world =
      penetration_mode.priorPenetrationDirection();

  // Cache
  PenetrationObjectCallbackCache<S> cache1;
  PenetrationObjectCallbackCache<S> cache2;
  GJKGeometryData<S> geom1;
  GJKGeometryData<S> geom2;
  Transform3<S> tf_geom_1;
  Transform3<S> tf_geom_2;

  // Make mpr
  const int mpr_max_iterations = 128;
  MPR<S> mpr(mpr_max_iterations, request.distanceTolerance());

  // Result processing functor
  auto process_result = [&](const Contact<S>& c, bool& keep_this,
                            bool& user_stop) -> void {
    // Get the object
    penetrationDistanceGetContactGJK(c, o1, tf1, true, geom1, tf_geom_1,
                                     cache1);
    penetrationDistanceGetContactGJK(c, o2, tf2, true, geom2, tf_geom_2,
                                     cache2);

    // Construct the new contact
    Contact<S> new_contact = c;
    bool is_reversed_o1o2 = (c.o1 == c.o2) ? false : (o1 == c.o2);
    if (is_reversed_o1o2) {
      penetrationDepthSwapContactObjectPair(new_contact);
    }

    // Invoke PD
    const bool is_incremental_penetration =
        (penetration_mode.penetration_type() ==
         CollisionPenetrationType::IncrementalMinimumPenetration);
    computePenetrationMPR(geom1, tf_geom_1, geom2, tf_geom_2,
                          shape2_escape_direction_world, mpr,
                          is_incremental_penetration, new_contact.pos,
                          new_contact.normal, new_contact.penetration_depth);
    result.addContact(new_contact);

    // Assign the flag
    user_stop = result.terminationConditionSatisfied(request.maxNumContacts());
    keep_this = false;
  };

  CollisionRequest<S> dummy_request = request;
  dummy_request.disablePenetration();
  CollisionResult<S> dummy_result(std::move(process_result));
  ::fcl::detail::collide(o1, tf1, o2, tf2, &solver, dummy_request,
                         dummy_result);
  return result.numContacts();
}

}  // namespace detail
}  // namespace fcl