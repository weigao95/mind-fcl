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
void randomOctreePairCollision(std::uint16_t bottom_half_shape = 2,
                               std::size_t test_n_points = 10,
                               std::size_t test_n_collision = 10,
                               bool compare_with_naive = true) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  auto octree1 = test::makeRandomPointsAsOctrees(
      scalar_resolution, bottom_half_shape, test_n_points);
  auto octree2 = test::makeRandomPointsAsOctrees(
      scalar_resolution, bottom_half_shape, test_n_points);

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  // Try solve
  const auto& raw_tree1 = octree1->raw_octree();
  const auto& raw_tree2 = octree2->raw_octree();
  detail::GJKSolver<S> gjk_solver;
  detail::CollisionSolverOctree2<S> solver(&gjk_solver);
  Transform3<S> tree1_pose, tree2_pose;
  CollisionRequest<S> request;
  request.useDefaultPenetration();
  request.setMaxContactCount(100000);
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, tree1_pose);
    test::generateRandomTransform(extent, tree2_pose);

    CollisionResult<S> result_i;
    solver.OctreePairIntersect(octree1.get(), octree2.get(), tree1_pose,
                               tree2_pose, request, result_i);

    // Verify contacts
    const auto& contacts = result_i.getContacts();
    for (const Contact<S>& contact_j : contacts) {
      EXPECT_TRUE(raw_tree1->isPointOccupied(contact_j.o1_bv.center()));
      EXPECT_TRUE(raw_tree2->isPointOccupied(contact_j.o2_bv.center()));

      Box<S> box1_local, box2_local;
      Transform3<S> box1_tf_local, box2_tf_local;
      constructBox(contact_j.o1_bv, tree1_pose, box1_local, box1_tf_local);
      constructBox(contact_j.o2_bv, tree2_pose, box2_local, box2_tf_local);
      CollisionResult<S> box_pair_result;
      fcl::collide(&box1_local, box1_tf_local, &box2_local, box2_tf_local,
                   request, box_pair_result);
      EXPECT_TRUE(box_pair_result.isCollision());
    }

    // Maybe too large
    if (!compare_with_naive) continue;

    // As leaf testing
    std::vector<Contact<S>> traverse_colliding_contacts;
    auto visit_fn = [&](const AABB<S>& bv) -> bool {
      Box<S> box1_local;
      Transform3<S> box1_tf_local;
      constructBox(bv, tree1_pose, box1_local, box1_tf_local);
      CollisionResult<S> box_octree_result;
      solver.ShapeOctreeIntersect(box1_local, octree2.get(), box1_tf_local,
                                  tree2_pose, request, box_octree_result);
      const auto& box_octree_contacts = box_octree_result.getContacts();
      for (const auto& new_contact : box_octree_contacts) {
        traverse_colliding_contacts.push_back(new_contact);
      }

      // Done
      return false;
    };
    octree1->visitLeafNodes(visit_fn);
    EXPECT_EQ(traverse_colliding_contacts.size(), contacts.size());
    // std::cout << "Random octree pair contact size " << contacts.size() <<
    // std::endl;
  }
}

}  // namespace fcl

GTEST_TEST(OctreePairCollision, RandomOctreeTest) {
  fcl::randomOctreePairCollision<float>(2, 50, 10, true);
  fcl::randomOctreePairCollision<double>(2, 50, 10, true);
  fcl::randomOctreePairCollision<float>(8, 4000, 10, true);
  fcl::randomOctreePairCollision<double>(8, 4000, 10, true);
  fcl::randomOctreePairCollision<float>(16, 10000, 10, false);
  fcl::randomOctreePairCollision<double>(16, 10000, 10, false);
  fcl::randomOctreePairCollision<float>(32, 100000, 10, false);
  fcl::randomOctreePairCollision<double>(32, 100000, 10, false);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}