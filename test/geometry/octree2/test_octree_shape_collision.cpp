//
// Created by mech-mind_gw on 3/26/2024.
//

#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void plane_xOy_OctreeWithBox(std::uint16_t bottom_half_shape = 2,
                             std::size_t test_n_collision = 10) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  auto tree = test::makePlane_xOy_AsOctree2(bottom_half_shape, scalar_resolution);
  Box<S> box(bottom_half_size * 0.3, bottom_half_size * 0.15,
             bottom_half_size * 0.2);

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  const auto& raw_tree = tree->raw_octree();
  detail::GJKSolver<S> gjk_solver;
  detail::CollisionSolverOctree2<S> solver(&gjk_solver);
  Transform3<S> octree_pose, box_pose;
  CollisionRequest<S> request;
  request.setMaxContactCount(100000);
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, octree_pose);
    test::generateRandomTransform(extent, box_pose);

    CollisionResult<S> result_i;
    solver.OctreeShapeIntersect(tree.get(), box, octree_pose, box_pose,
                                request, result_i);

    // Verify contacts
    const auto& contacts = result_i.getContacts();
    for (const Contact<S>& contact_j : contacts) {
      const Vector3<S> voxel_center = contact_j.o1_bv.center();
      EXPECT_TRUE(raw_tree->isPointOccupied(voxel_center));

      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(contact_j.o1_bv, octree_pose, box_local, box_tf_local);
      CollisionResult<S> box_pair_result;
      fcl::collide(&box_local, box_tf_local, &box, box_pose, request,
                   box_pair_result);
      EXPECT_TRUE(box_pair_result.isCollision());
    }

    // Very by traverse
    std::vector<AABB<S>> traverse_colliding_AABB;
    auto visit_fn = [&](const AABB<S>& bv) -> bool {
      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(bv, octree_pose, box_local, box_tf_local);
      CollisionResult<S> box_pair_result;
      fcl::collide(&box_local, box_tf_local, &box, box_pose, request,
                   box_pair_result);
      if (box_pair_result.isCollision()) {
        traverse_colliding_AABB.push_back(bv);
      }

      // Continue
      return false;
    };
    tree->visitLeafNodes(visit_fn);
    EXPECT_EQ(contacts.size(), traverse_colliding_AABB.size());
  }
}

template <typename S>
void randomOctreeCollisionWithBox(std::uint16_t bottom_half_shape = 2,
                                  std::size_t test_n_points = 10,
                                  std::size_t test_n_collision = 10) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  std::shared_ptr<const fcl::Octree2CollisionGeometry<S>> octree = nullptr;
  {
    using namespace octree2;
    Vector3<S> resolution(scalar_resolution, scalar_resolution,
                                 scalar_resolution);
    std::vector<Vector3<S>> points_inserted;
    for (std::size_t i = 0; i < test_n_points; i++) {
      Vector3<S> point_i;
      point_i.setRandom();
      point_i *= (0.99 * bottom_half_size);
      points_inserted.push_back(point_i);
    }

    // Insert into tree
    auto point_fn = [&points_inserted](int index, S& x, S& y, S& z) -> void {
      const auto point = points_inserted[index];
      x = point.x();
      y = point.y();
      z = point.z();
    };
    auto tree = std::make_shared<Octree<S>>(resolution, bottom_half_shape);
    tree->rebuildTree(point_fn, points_inserted.size());

    auto mutable_octree = std::make_shared<Octree2CollisionGeometry<S>>(tree);
    mutable_octree->computeLocalAABB();
    octree = mutable_octree;
  }

  // Make a box
  Box<S> box(bottom_half_size * 0.3, bottom_half_size * 0.15,
             bottom_half_size * 0.2);

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  // Try solve
  const auto& raw_tree = octree->raw_octree();
  detail::GJKSolver<S> gjk_solver;
  detail::CollisionSolverOctree2<S> solver(&gjk_solver);
  Transform3<S> octree_pose, box_pose;
  CollisionRequest<S> request;
  request.setMaxContactCount(100000);
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, octree_pose);
    test::generateRandomTransform(extent, box_pose);

    CollisionResult<S> result_i;
    solver.OctreeShapeIntersect(octree.get(), box, octree_pose, box_pose,
                                request, result_i);

    // Verify contacts
    const auto& contacts = result_i.getContacts();
    for (const Contact<S>& contact_j : contacts) {
      const Vector3<S> voxel_center = contact_j.o1_bv.center();
      EXPECT_TRUE(raw_tree->isPointOccupied(voxel_center));

      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(contact_j.o1_bv, octree_pose, box_local, box_tf_local);
      CollisionResult<S> box_pair_result;
      fcl::collide(&box_local, box_tf_local, &box, box_pose, request,
                   box_pair_result);
      EXPECT_TRUE(box_pair_result.isCollision());
    }

    // Very by traverse
    std::vector<AABB<S>> traverse_colliding_AABB;
    auto visit_fn = [&](const AABB<S>& bv) -> bool {
      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(bv, octree_pose, box_local, box_tf_local);
      CollisionResult<S> box_pair_result;
      fcl::collide(&box_local, box_tf_local, &box, box_pose, request,
                   box_pair_result);
      if (box_pair_result.isCollision()) {
        traverse_colliding_AABB.push_back(bv);
      }

      // Continue
      return false;
    };
    octree->visitLeafNodes(visit_fn);
    EXPECT_EQ(contacts.size(), traverse_colliding_AABB.size());
  }
}

}  // namespace fcl

GTEST_TEST(Octree2ShapeCollision, RandomBoxTest) {
  fcl::randomOctreeCollisionWithBox<float>(2, 20, 20);
  fcl::randomOctreeCollisionWithBox<double>(2, 20, 20);
  fcl::randomOctreeCollisionWithBox<float>(64, 1000 * 100, 10);
  fcl::randomOctreeCollisionWithBox<double>(64, 1000 * 100, 10);
}

GTEST_TEST(Octree2ShapeCollision, Plane_xOy_BoxTest) {
  fcl::plane_xOy_OctreeWithBox<float>(2, 20);
  fcl::plane_xOy_OctreeWithBox<double>(2, 20);
  fcl::plane_xOy_OctreeWithBox<float>(1024, 20);
  fcl::plane_xOy_OctreeWithBox<double>(1024, 20);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}