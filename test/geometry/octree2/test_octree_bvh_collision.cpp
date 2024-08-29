//
// Created by mech-mind_gw on 3/27/2024.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_collision_geometry.h"
#include "fcl/narrowphase/detail/traversal/octree2/octree2_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void randomOctreeBoxMeshCollision(std::uint16_t bottom_half_shape = 2,
                                  std::size_t test_n_points = 10,
                                  std::size_t test_n_collision = 10,
                                  bool compare_with_naive = true) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  auto octree1 = test::makeRandomPointsAsOctrees(
      scalar_resolution, bottom_half_shape, test_n_points);

  // Make a box
  Box<S> box(0.4 * bottom_half_size, 0.4 * bottom_half_size,
             0.4 * bottom_half_size);
  auto mesh = test::generateBoxBVHModel<OBB<S>>(box);

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  // Try solve
  const auto& raw_tree1 = octree1->raw_octree();
  detail::GJKSolver<S> gjk_solver;
  detail::CollisionSolverOctree2<S> solver(&gjk_solver);
  Transform3<S> tree1_pose, mesh2_pose;
  CollisionRequest<S> request;
  request.setMaxContactCount(100000);
  for (std::size_t i = 0; i < test_n_collision; i++) {
    test::generateRandomTransform(extent, tree1_pose);
    test::generateRandomTransform(extent, mesh2_pose);

    CollisionResult<S> result_i;
    solver.OctreeBVHIntersect(octree1.get(), mesh.get(), tree1_pose, mesh2_pose,
                              request, result_i);

    // Verify contacts
    const auto& contacts = result_i.getContacts();
    for (const Contact<S>& contact_j : contacts) {
      EXPECT_TRUE(raw_tree1->isPointOccupied(contact_j.o1_bv.center()));

      // Make the box
      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(contact_j.o1_bv, tree1_pose, box_local, box_tf_local);

      // Make the triangle
      auto bvh_primitive_id = contact_j.b2;
      const MeshSimplex& tri_id = mesh->simplex_indices()[bvh_primitive_id];
      const Vector3<S>& p1 = mesh->getVertex(tri_id[0]);
      const Vector3<S>& p2 = mesh->getVertex(tri_id[1]);
      const Vector3<S>& p3 = mesh->getVertex(tri_id[2]);
      bool intersect = gjk_solver.shapeTriangleIntersect(
          box_local, box_tf_local, p1, p2, p3, mesh2_pose);
      EXPECT_TRUE(intersect);
    }

    if (!compare_with_naive) continue;

    // Collision checking with naive
    std::size_t naive_count = 0;
    CollisionResult<S> result1;
    for (int primitive_id = 0; primitive_id < mesh->num_simplex();
         primitive_id++) {
      const MeshSimplex& tri_id = mesh->simplex_indices()[primitive_id];
      const Vector3<S>& p1 = mesh->getVertex(tri_id[0]);
      const Vector3<S>& p2 = mesh->getVertex(tri_id[1]);
      const Vector3<S>& p3 = mesh->getVertex(tri_id[2]);
      fcl::TriangleP<S> triangle_p(p1, p2, p3);
      result1.clear();
      solver.ShapeOctreeIntersect(triangle_p, octree1.get(), mesh2_pose,
                                  tree1_pose, request, result1);
      naive_count += result1.numContacts();
    }

    EXPECT_EQ(naive_count, contacts.size());
  }
}

}  // namespace fcl

GTEST_TEST(OctreeMeshCollision, RandomOctreeBoxMesh) {
  fcl::randomOctreeBoxMeshCollision<float>(2, 50, 10, true);
  fcl::randomOctreeBoxMeshCollision<double>(2, 50, 10, true);
  fcl::randomOctreeBoxMeshCollision<float>(8, 4000, 10, true);
  fcl::randomOctreeBoxMeshCollision<double>(8, 4000, 10, true);
  fcl::randomOctreeBoxMeshCollision<float>(16, 10000, 10, true);
  fcl::randomOctreeBoxMeshCollision<double>(16, 10000, 10, true);
  fcl::randomOctreeBoxMeshCollision<float>(32, 100000, 10, true);
  fcl::randomOctreeBoxMeshCollision<double>(32, 100000, 10, true);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}