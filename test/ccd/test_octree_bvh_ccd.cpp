//
// Created by Wei Gao on 2024/6/17.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void randomOctreeBoxMeshCCD(std::uint16_t bottom_half_shape = 2,
                            std::size_t test_n_points = 10,
                            std::size_t test_n_collision = 10,
                            bool compare_with_naive = true) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  auto octree1 = test::makeRandomPointsAsOctrees(
      scalar_resolution, bottom_half_shape, test_n_points);

  // Make a box mesh
  std::shared_ptr<const BVHModel<OBB<S>>> mesh;
  {
    Box<S> box(0.4 * bottom_half_size, 0.4 * bottom_half_size,
               0.4 * bottom_half_size);
    mesh = test::generateBoxBVHModel<OBB<S>>(box);
  }

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  // Try solve
  const auto& raw_tree1 = octree1->raw_octree();
  detail::TranslationalDisplacementOctreeSolver<S> solver;
  Transform3<S> tree1_pose, mesh2_pose;
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 100000;
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, tree1_pose);
    test::generateRandomTransform(extent, mesh2_pose);

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
    solver.RunOctreeObbBVH(octree1.get(), tree1_pose, tree1_displacement,
                           mesh.get(), mesh2_pose, request, result_i);

    // Verify contacts
    const auto& contacts = result_i.raw_contacts();
    for (const ContinuousCollisionContact<S>& contact_j : contacts) {
      EXPECT_TRUE(raw_tree1->isPointOccupied(contact_j.o1_bv.center()));

      // Make the box
      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(contact_j.o1_bv, tree1_pose, box_local, box_tf_local);
      box_local.computeLocalAABB();

      // Make the triangle
      auto bvh_primitive_id = contact_j.b2;
      const MeshSimplex& tri_id = mesh->simplex_indices()[bvh_primitive_id];
      const Vector3<S>& p1 = mesh->getVertex(tri_id[0]);
      const Vector3<S>& p2 = mesh->getVertex(tri_id[1]);
      const Vector3<S>& p3 = mesh->getVertex(tri_id[2]);
      Simplex<S> triangle(p1, p2, p3);
      ContinuousCollisionResult<S> box_result;
      detail::ContinuousContactMeta<S> box_contact_meta;
      using ShapePairSolver = detail::ShapePairTranslationalCollisionSolver<S>;
      ShapePairSolver::RunShapeSimplex(box_local, box_tf_local,
                                       tree1_displacement, triangle, mesh2_pose,
                                       box_contact_meta, request, box_result);
      EXPECT_TRUE(box_result.num_contacts() > 0);
    }

    if (!compare_with_naive) continue;

    // Collision checking with naive
    std::size_t naive_count = 0;
    // std::size_t naive_count2 = 0;
    for (int primitive_id = 0; primitive_id < mesh->num_simplex();
         primitive_id++) {
      const MeshSimplex& tri_id = mesh->simplex_indices()[primitive_id];
      const Vector3<S>& p1 = mesh->getVertex(tri_id[0]);
      const Vector3<S>& p2 = mesh->getVertex(tri_id[1]);
      const Vector3<S>& p3 = mesh->getVertex(tri_id[2]);
      fcl::Simplex<S> simplex_primitive(p1, p2, p3);

      // By visitor
      std::vector<ContinuousCollisionContact<S>> contact_for_simplex;
      std::size_t count_for_this_primitive = 0;
      auto visitor = [&](const AABB<S>& aabb) -> bool {
        // Make the box
        Box<S> box_local;
        Transform3<S> box_tf_local;
        constructBox(aabb, tree1_pose, box_local, box_tf_local);
        box_local.computeLocalAABB();

        ContinuousCollisionResult<S> box_result;
        detail::ContinuousContactMeta<S> box_contact_meta;
        box_contact_meta.o1_bv = aabb;
        using ShapePairSolver =
            detail::ShapePairTranslationalCollisionSolver<S>;
        ShapePairSolver::RunShapeSimplex(
            box_local, box_tf_local, tree1_displacement, simplex_primitive,
            mesh2_pose, box_contact_meta, request, box_result);

        count_for_this_primitive += box_result.num_contacts();
        for (const auto& elem : box_result.raw_contacts()) {
          contact_for_simplex.push_back(elem);
        }

        // Not finish
        return false;
      };
      octree1->visitLeafNodes(visitor);

      // Update naive_count
      naive_count += count_for_this_primitive;
    }

    EXPECT_EQ(naive_count, contacts.size());
  }
}

}  // namespace fcl

GTEST_TEST(OctreeMeshCollisionCCD, RandomOctreeBoxMesh) {
  // fcl::testInstancePair<float>();
  fcl::randomOctreeBoxMeshCCD<float>(2, 50, 10, true);
  fcl::randomOctreeBoxMeshCCD<double>(2, 50, 10, true);
  fcl::randomOctreeBoxMeshCCD<float>(8, 4000, 10, true);
  fcl::randomOctreeBoxMeshCCD<double>(8, 4000, 10, true);
  fcl::randomOctreeBoxMeshCCD<float>(16, 10000, 10, true);
  fcl::randomOctreeBoxMeshCCD<double>(16, 10000, 10, true);
  fcl::randomOctreeBoxMeshCCD<float>(32, 100000, 10, true);
  fcl::randomOctreeBoxMeshCCD<double>(32, 100000, 10, true);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}