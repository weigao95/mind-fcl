//
// Created by wei on 24-6-16.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/continuous_collision.h"
#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void randomOctreePairCollisionCCD(std::uint16_t bottom_half_shape = 2,
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
  detail::TranslationalDisplacementOctreeSolver<S> solver;
  Transform3<S> tree1_pose, tree2_pose;
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 100000;
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, tree1_pose);
    test::generateRandomTransform(extent, tree2_pose);

    // Random displacement
    TranslationalDisplacement<S> tree1_displacement;
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);
      Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

      // Setup the axis
      tree1_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
      tree1_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    ContinuousCollisionResult<S> result_i;
    solver.RunOctreePair(octree1.get(), tree1_pose, tree1_displacement,
                         octree2.get(), tree2_pose, request, result_i);

    // Verify contacts
    const auto& contacts = result_i.raw_contacts();
    for (const ContinuousCollisionContact<S>& contact_j : contacts) {
      EXPECT_TRUE(raw_tree1->isPointOccupied(contact_j.o1_bv.center()));
      EXPECT_TRUE(raw_tree2->isPointOccupied(contact_j.o2_bv.center()));

      Box<S> box1_local, box2_local;
      Transform3<S> box1_tf_local, box2_tf_local;
      constructBox(contact_j.o1_bv, tree1_pose, box1_local, box1_tf_local);
      constructBox(contact_j.o2_bv, tree2_pose, box2_local, box2_tf_local);
      box1_local.computeLocalAABB();
      box2_local.computeLocalAABB();
      ContinuousCollisionResult<S> box_pair_result;
      fcl::translational_ccd<S>(&box1_local, box1_tf_local, tree1_displacement,
                                &box2_local, box2_tf_local, request,
                                box_pair_result);
      /*std::cout << "Box1 pose" << std::endl;
      std::cout << box1_tf_local.matrix() << std::endl;
      std::cout << "Box1 size" << std::endl;
      std::cout << box1_local.side << std::endl;
      std::cout << "Box2 pose" << std::endl;
      std::cout << box2_tf_local.matrix() << std::endl;
      std::cout << "Box2 size" << std::endl;
      std::cout << box2_local.side << std::endl;*/
      if (box_pair_result.num_contacts() == 0) {
        ContinuousCollisionResult<S> tmp_result_i;
        solver.RunOctreePair(octree1.get(), tree1_pose, tree1_displacement,
                             octree2.get(), tree2_pose, request, tmp_result_i);
      }
      EXPECT_TRUE(box_pair_result.num_contacts() > 0);
    }

    // Maybe too large
    if (!compare_with_naive) continue;

    // As leaf testing
    std::vector<ContinuousCollisionContact<S>> traverse_colliding_contacts;
    auto visit_fn = [&](const AABB<S>& bv) -> bool {
      Box<S> box1_local;
      Transform3<S> box1_tf_local;
      constructBox(bv, tree1_pose, box1_local, box1_tf_local);
      box1_local.computeLocalAABB();
      ContinuousCollisionResult<S> box_octree_result;
      solver.RunShapeOctree(&box1_local, box1_tf_local, tree1_displacement,
                            octree2.get(), tree2_pose, request,
                            box_octree_result);
      const auto& box_octree_contacts = box_octree_result.raw_contacts();
      for (const auto& new_contact : box_octree_contacts) {
        traverse_colliding_contacts.push_back(new_contact);
      }

      // Done
      return false;
    };
    octree1->visitLeafNodes(visit_fn);
    EXPECT_EQ(traverse_colliding_contacts.size(), contacts.size());
    // std::cout << "Random translational octree pair contact size "
    //           << contacts.size() << std::endl;
  }
}

}  // namespace fcl

GTEST_TEST(OctreePairTranslationalCollision, RandomOctreeTest) {
  fcl::randomOctreePairCollisionCCD<float>(2, 50, 10, true);
  fcl::randomOctreePairCollisionCCD<double>(2, 50, 10, true);
  fcl::randomOctreePairCollisionCCD<float>(8, 4000, 10, true);
  fcl::randomOctreePairCollisionCCD<double>(8, 4000, 10, true);
  fcl::randomOctreePairCollisionCCD<float>(16, 10000, 10, false);
  fcl::randomOctreePairCollisionCCD<double>(16, 10000, 10, false);
  fcl::randomOctreePairCollisionCCD<float>(32, 100000, 10, false);
  fcl::randomOctreePairCollisionCCD<double>(32, 100000, 10, false);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}