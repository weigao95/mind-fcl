//
// Created by Wei Gao on 2024/6/17.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/narrowphase/detail/ccd/octree2_ccd_solver.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void testInstancePair() {
  Box<S> box(0.400000006, 0.400000006, 0.400000006);
  const auto val = 0.159999996;
  auto a = Vector3<S>(val, val, val);
  auto b = Vector3<S>(val, -val, -val);
  auto c = Vector3<S>(val, val, -val);
  TriangleP<S> triangle(a, b, c);
  Transform3<S> box_tf, tri_tf;
  box_tf.setIdentity();
  box_tf.translation() = Vector3<S>(0.2, -0.2, 0.2);

  tri_tf.setIdentity();
  // tri_tf.translation() = Vector3<S>(0.000, 0.000, 0.001);

  const auto scalar_disp = 0.1;
  TranslationalDisplacement<S> box_displacement;
  box_displacement.unit_axis_in_shape1 = Vector3<S>::UnitZ();
  box_displacement.scalar_displacement = scalar_disp;
  TranslationalDisplacement<S> equiv_tri_displacement;
  equiv_tri_displacement.unit_axis_in_shape1 = -Vector3<S>::UnitZ();
  equiv_tri_displacement.scalar_displacement = scalar_disp;

  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 100000;
  using ShapePairSolver = detail::ShapePairTranslationalCollisionSolver<S>;
  ContinuousCollisionResult<S> box_result, tri_result;
  ShapePairSolver::RunShapePair(&box, box_tf, box_displacement, &triangle,
                                tri_tf, request, box_result);
  ShapePairSolver::RunShapePair(&triangle, tri_tf, equiv_tri_displacement, &box,
                                box_tf, request, tri_result);
  EXPECT_EQ(box_result.num_contacts(), tri_result.num_contacts());

  Simplex<S> simplex(a, b, c);
  detail::ContinuousContactMeta<S> meta;
  ContinuousCollisionResult<S> simplex_result;
  ShapePairSolver::RunShapeSimplex(box, box_tf, box_displacement, simplex,
                                   tri_tf, meta, request, simplex_result);
  EXPECT_EQ(box_result.num_contacts(), simplex_result.num_contacts());
}

template <typename S>
void testTrianglesInMesh(std::uint16_t bottom_half_shape = 2,
                         std::size_t test_n_points = 10,
                         std::size_t test_n_collision = 10) {
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
  detail::TranslationalDisplacementOctreeSolver<S> solver;
  Transform3<S> tree1_pose, mesh2_pose;
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 100000000;
  for (std::size_t test_i = 0; test_i < test_n_collision; test_i++) {
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

    TranslationalDisplacement<S> mesh_displacement;
    mesh_displacement.scalar_displacement =
        tree1_displacement.scalar_displacement;
    mesh_displacement.unit_axis_in_shape1 =
        mesh2_pose.linear().transpose() *
        (tree1_pose.linear() * (-tree1_displacement.unit_axis_in_shape1));

    // std::size_t naive_count2 = 0;
    for (int primitive_id = 0; primitive_id < mesh->num_simplex();
         primitive_id++) {
      const MeshSimplex& tri_id = mesh->simplex_indices()[primitive_id];
      const Vector3<S>& p1 = mesh->getVertex(tri_id[0]);
      const Vector3<S>& p2 = mesh->getVertex(tri_id[1]);
      const Vector3<S>& p3 = mesh->getVertex(tri_id[2]);

      // As triangle
      ContinuousCollisionResult<S> octree_triangle_result;
      fcl::TriangleP<S> triangle_p(p1, p2, p3);
      triangle_p.computeLocalAABB();
      // std::cout << "Octree Shape " << primitive_id << std::endl;
      solver.template RunOctreeShape<TriangleP<S>>(
          octree1.get(), tree1_pose, tree1_displacement, &triangle_p,
          mesh2_pose, request, octree_triangle_result);

      // std::cout << "Shape Octree " << primitive_id << std::endl;
      ContinuousCollisionResult<S> triangle_octree_result;
      solver.template RunShapeOctree<TriangleP<S>>(
          &triangle_p, mesh2_pose, mesh_displacement, octree1.get(), tree1_pose,
          request, triangle_octree_result);
      EXPECT_EQ(octree_triangle_result.num_contacts(),
                triangle_octree_result.num_contacts());

      // As simplex
      fcl::Simplex<S> simplex_primitive(p1, p2, p3);
      std::vector<ContinuousCollisionContact<S>> contact_for_simplex;
      std::size_t count_for_this_primitive = 0;
      auto visitor = [&](const AABB<S>& aabb) -> bool {
        // Make the box
        Box<S> box_local;
        Transform3<S> box_tf_local;
        constructBox(aabb, tree1_pose, box_local, box_tf_local);
        box.computeLocalAABB();

        ContinuousCollisionResult<S> box_result;
        EXPECT_EQ(box_result.num_contacts(), 0u);
        detail::ContinuousContactMeta<S> box_contact_meta;
        box_contact_meta.o1_bv = aabb;
        using ShapePairSolver =
            detail::ShapePairTranslationalCollisionSolver<S>;
        ShapePairSolver::RunShapePair(
            box_local, box_tf_local, tree1_displacement, triangle_p, mesh2_pose,
            box_contact_meta, request, box_result);

        count_for_this_primitive += box_result.num_contacts();
        for (const auto& elem : box_result.raw_contacts()) {
          contact_for_simplex.push_back(elem);
        }

        // Not finish
        return false;
      };
      octree1->visitLeafNodes(visitor);

      EXPECT_EQ(count_for_this_primitive,
                octree_triangle_result.num_contacts());
    }
  }
}

template <typename S>
void testRandomOctreeWithBox(std::uint16_t bottom_half_shape = 2,
                             std::size_t test_n_points = 10,
                             std::size_t test_n_collision = 10) {
  // Make the tree
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  auto octree1 = test::makeRandomPointsAsOctrees(
      scalar_resolution, bottom_half_shape, test_n_points);
  Box<S> box(bottom_half_size * 0.3, bottom_half_size * 0.15,
             bottom_half_size * 0.2);
  box.computeLocalAABB();

  // Random pose
  const S extent_scalar = bottom_half_size * 0.3;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};

  // Try solve
  detail::TranslationalDisplacementOctreeSolver<S> solver;
  Transform3<S> octree_pose, box_pose;
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 100000000;
  for (std::size_t test_i = 0; test_i < test_n_collision; test_i++) {
    test::generateRandomTransform(extent, octree_pose);
    test::generateRandomTransform(extent, box_pose);

    // Random displacement
    TranslationalDisplacement<S> octree_displacement;
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);
      Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

      // Setup the axis
      octree_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
      octree_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    TranslationalDisplacement<S> box_equiv_displacement;
    box_equiv_displacement.scalar_displacement =
        octree_displacement.scalar_displacement;
    box_equiv_displacement.unit_axis_in_shape1 =
        box_pose.linear().transpose() *
        (octree_pose.linear() * (-octree_displacement.unit_axis_in_shape1));

    ContinuousCollisionResult<S> octree_shape_result;
    solver.template RunOctreeShape<Box<S>>(octree1.get(), octree_pose,
                                           octree_displacement, &box, box_pose,
                                           request, octree_shape_result);
    ContinuousCollisionResult<S> shape_octree_result;
    solver.template RunShapeOctree<Box<S>>(
        &box, box_pose, box_equiv_displacement, octree1.get(), octree_pose,
        request, shape_octree_result);
    EXPECT_EQ(octree_shape_result.num_contacts(),
              shape_octree_result.num_contacts());

    std::size_t count_naive = 0;
    auto visitor = [&](const AABB<S>& aabb) -> bool {
      // Make the box
      Box<S> box_local;
      Transform3<S> box_tf_local;
      constructBox(aabb, octree_pose, box_local, box_tf_local);
      box_local.computeLocalAABB();

      ContinuousCollisionResult<S> box_result;
      EXPECT_EQ(box_result.num_contacts(), 0u);
      detail::ContinuousContactMeta<S> box_contact_meta;
      box_contact_meta.o1_bv = aabb;
      using ShapePairSolver =
          detail::ShapePairTranslationalCollisionSolver<S>;
      ShapePairSolver::RunShapePair(
          box_local, box_tf_local, octree_displacement, box, box_pose,
          box_contact_meta, request, box_result);
      count_naive += box_result.num_contacts();

      // Not finish
      return false;
    };
    octree1->visitLeafNodes(visitor);
    EXPECT_EQ(count_naive, octree_shape_result.num_contacts());
  }
}

}  // namespace fcl

GTEST_TEST(OctreeShapeCollisionCCD, SingularCase) {
  // fcl::testInstancePair<float>();
}

GTEST_TEST(OctreeShapeCollisionCCD, TrianglesInMesh) {
  fcl::testTrianglesInMesh<float>(2, 50, 100);
  fcl::testTrianglesInMesh<double>(2, 50, 100);
  fcl::testTrianglesInMesh<float>(8, 4000, 10);
  fcl::testTrianglesInMesh<double>(8, 4000, 10);
}

GTEST_TEST(OctreeShapeCollisionCCD, RandomOctreeBox) {
  fcl::testRandomOctreeWithBox<float>(2, 50, 100);
  fcl::testRandomOctreeWithBox<double>(2, 50, 100);
  fcl::testRandomOctreeWithBox<float>(8, 4000, 10);
  fcl::testRandomOctreeWithBox<double>(8, 4000, 10);
}

int main(int argc, char* argv[]) {
  srand(1);
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}