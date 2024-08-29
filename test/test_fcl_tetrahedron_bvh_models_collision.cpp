#include <gtest/gtest.h>

#include <iostream>

#include "eigen_matrix_compare.h"
#include "make_box_mesh.h"
#include "make_cylinder_mesh.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace test {

GTEST_TEST(MakeBoxBvhMeshTest, CalcSequentialIndex) {
  const Vector3<int> num_vertices(3, 2, 5);
  EXPECT_EQ(28, CalcSequentialIndex(2, 1, 3, num_vertices));
}

GTEST_TEST(MakeBoxBvhMeshTest, GenerateVertices) {
  // Set up a box [-1,1]x[-2,2]x[-3,3] whose corners have integer coordinates.
  const Box<double> box(2.0, 4.0, 6.0);
  // Request the number of vertices so that the vertices have integer
  // coordinates {-1,0,1} x {-2,-1,0,1,2} x {-3,-2,-1,0,1,2,3}.
  // Vertex[i][j][k] should have coordinates (i-1, j-2, k-3) for
  // 0 ≤ i < 3, 0 ≤ j < 5, 0 ≤ k < 7.
  const Vector3<int> num_vertices{3, 5, 7};

  auto vertices = GenerateVertices<double>(box, num_vertices);

  EXPECT_EQ(105, static_cast<int>(vertices.size()));
  for (int i = 0; i < num_vertices.x(); ++i) {
    for (int j = 0; j < num_vertices.y(); ++j) {
      for (int k = 0; k < num_vertices.z(); ++k) {
        int sequential_index = CalcSequentialIndex(i, j, k, num_vertices);
        Vector3<double> expect_r_MV = Vector3<double>(i - 1, j - 2, k - 3);
        Vector3<double> r_MV = vertices[sequential_index];
        EXPECT_TRUE(CompareMatrices(expect_r_MV, r_MV))
            << "Incorrect vertex position.";
      }
    }
  }
}

GTEST_TEST(MakeBoxBvhMeshTest, AddSixTetrahedraOfCell) {
  const Vector3<int> lowest(1, 2, 3);
  const Vector3<int> num_vertices(3, 4, 5);
  std::vector<MeshSimplex> elements;

  AddSixTetrahedraOfCell(lowest, num_vertices, &elements);
  EXPECT_EQ(6, static_cast<int>(elements.size()));

  // In a 3x4x5 grid of vertices, the vertex with (i,j,k)-index = (1,2,3) has
  // its sequential index 33. This picture shows how the rectangular cell
  // with its lowest vertex v₃₃ looks like.
  //
  //               v₃₄     v₃₉
  //               ●------●
  //              /|     /|
  //             / | v₅₉/ |
  //        v₅₄ ●------●  |
  //            |  |   |  |
  //            |  ●---|--● v₃₈
  //            | /v₃₃ | /
  //            |/     |/
  //    +K  v₅₃ ●------● v₅₈
  //     |
  //     |
  //     o------+J
  //    /
  //   /
  // +I
  //
  // This table has the expected six tetrahedra of the rectangular cell.
  // They share the main diagonal v₃₃v₅₉.
  const int expect_elements[6][4]{// clang-format off
                                  {33, 59, 53, 58},
                                  {33, 59, 58, 38},
                                  {33, 59, 38, 39},
                                  {33, 59, 39, 34},
                                  {33, 59, 34, 54},
                                  {33, 59, 54, 53}};
  // clang-format on
  for (int e = 0; e < 6; ++e)
    for (int v = 0; v < 4; ++v)
      EXPECT_EQ(expect_elements[e][v], elements[e][v]);
}

GTEST_TEST(MakeBoxBvhMeshTest, GenerateElements) {
  const Vector3<int> num_vertices{3, 5, 7};
  const int expect_total_num_vertex =
      num_vertices.x() * num_vertices.y() * num_vertices.z();

  const Vector3<int> num_cell = num_vertices - Vector3<int>::Ones();
  const int expect_num_cell = num_cell.x() * num_cell.y() * num_cell.z();
  const int expect_num_element = 6 * expect_num_cell;

  auto elements = GenerateElements(num_vertices);

  EXPECT_EQ(expect_num_element, static_cast<int>(elements.size()));
  // TODO(DamrongGuoy): Find a better way to test `elements`. Currently we
  //  only test that each tetrahedron uses vertices with indices in the range
  //  [0, expect_total_num_vertex). Perhaps check Euler characteristic,
  //  i.e., #vertex - #edge + #triangle - #tetrahedron = 1.
  for (const auto& simplex : elements) {
    for (int v = 0; v < 4; ++v) {
      EXPECT_GE(simplex[v], 0);
      EXPECT_LT(simplex[v], expect_total_num_vertex);
    }
  }
}

GTEST_TEST(MakeBoxBvhMeshTest, GenerateMesh) {
  const Box<double> box(0.2, 0.4, 0.8);
  auto box_mesh = MakeBoxBVHTetrahedronModel<fcl::OBB<double>>(box, 0.1);

  const int rectangular_cells = 2 * 4 * 8;
  const int tetrahedra_per_cell = 6;
  const int expect_num_tetrahedra = rectangular_cells * tetrahedra_per_cell;
  EXPECT_EQ(expect_num_tetrahedra, box_mesh.num_simplex());

  const int expect_num_vertices = 3 * 5 * 9;
  EXPECT_EQ(expect_num_vertices, box_mesh.num_vertices());
}

GTEST_TEST(MakeCylinderVolumeMesh, CoarsestMesh) {
  const double radius = 1.0;
  const double length = 2.0;
  // A `resolution_hint` greater than √2 times the radius of the cylinder
  // should give the coarsest mesh. We use a scaling factor larger than
  // √2 = 1.41421356... in the sixth decimal digit. It should give a mesh of
  // rectangular prism with 24 tetrahedra.
  const double resolution_hint_above = 1.414214 * radius;
  auto mesh_coarse = MakeCylinderBVHTetrahedronModel<fcl::OBB<double>>(
      Cylinder<double>(radius, length), resolution_hint_above);
  EXPECT_EQ(24, mesh_coarse.num_simplex());

  // A `resolution_hint` slightly below √2 times the radius of the cylinder
  // should give the next refined mesh.  We use a scaling factor smaller than
  // √2 = 1.41421356... in the sixth decimal digit. It should give a mesh of
  // rectangular prism with 8 * 24 = 192 tetrahedra.
  const double resolution_hint_below = 1.414213 * radius;
  auto mesh_fine = MakeCylinderBVHTetrahedronModel<fcl::OBB<double>>(
      Cylinder<double>(radius, length), resolution_hint_below);
  EXPECT_EQ(192, mesh_fine.num_simplex());
}

GTEST_TEST(MakeCylinderVolumeMesh, CoarsestMesh2) {
  const double height = 2;
  const double radius = 1;
  double resolution_hint = 2.0;  // Hint to coarsest mesh.
  auto mesh0 = MakeCylinderBVHTetrahedronModel<fcl::OBB<double>>(
      Cylinder<double>(radius, height), resolution_hint);

  for (int level = 1; level < 6; ++level) {
    resolution_hint /= 2.0;
    auto mesh = MakeCylinderBVHTetrahedronModel<fcl::OBB<double>>(
        Cylinder<double>(radius, height), resolution_hint);

    // Verify the correct size. There are initially 24 tetrahedra that each
    // split into 8 sub tetrahedra.
    const int num_tetrahedra = 24 * std::pow(8, level);
    EXPECT_EQ(mesh.num_simplex(), num_tetrahedra);
  }
}

template <typename S>
int isOverlap(const std::vector<Tetrahedron<S>>& tetrahedrons,
              const Transform3<S>& tf1, const Transform3<S>& tf2) {
  detail::GJKSolver<S> gjk_solver;
  int collision_count = 0;
  for (const auto& tet1 : tetrahedrons) {
    for (const auto& tet2 : tetrahedrons) {
      const bool is_overlap =
          gjk_solver.shapeIntersect(tet1, tf1, tet2, tf2, nullptr);
      if (is_overlap) {
        collision_count++;
      }
    }
  }
  return collision_count;
}

struct int_pair_hash {
  std::size_t operator()(const std::pair<int, int>& key) const {
    std::size_t a = std::hash<int>()(key.first);
    std::size_t b = std::hash<int>()(key.second);
    return a >= b ? a * a + a + b : a + b * b;
  }
};

template <typename S>
using CollisionResultMap =
    std::unordered_map<std::pair<int, int>, CollisionPenetrationContactData<S>,
                       int_pair_hash>;

template <typename S>
bool isOverlap(const std::vector<Tetrahedron<S>>& tetrahedrons,
               const Transform3<S>& tf1, const Transform3<S>& tf2,
               CollisionResultMap<S>& collision_result) {
  detail::GJKSolver<S> gjk_solver;
  collision_result.clear();
  for (int i = 0; i < static_cast<int>(tetrahedrons.size()); i++) {
    const auto& tet1 = tetrahedrons[i];
    for (int j = 0; j < static_cast<int>(tetrahedrons.size()); j++) {
      const auto& tet2 = tetrahedrons[j];
      CollisionPenetrationContactData<S> contacts;
      const bool is_overlap =
          gjk_solver.shapeIntersect(tet1, tf1, tet2, tf2, &contacts);
      if (is_overlap) {
        collision_result[std::make_pair(i, j)] = std::move(contacts);
      }
    }
  }
  return !collision_result.empty();
}

template <typename S>
bool isOverlap(const std::vector<Tetrahedron<S>>& tetrahedrons1,
               const std::vector<Tetrahedron<S>>& tetrahedrons2,
               const Transform3<S>& tf1, const Transform3<S>& tf2,
               CollisionResultMap<S>& collision_result) {
  detail::GJKSolver<S> gjk_solver;
  collision_result.clear();
  for (int i = 0; i < static_cast<int>(tetrahedrons1.size()); i++) {
    const auto& tet1 = tetrahedrons1[i];
    for (int j = 0; j < static_cast<int>(tetrahedrons2.size()); j++) {
      const auto& tet2 = tetrahedrons2[j];
      CollisionPenetrationContactData<S> contacts;
      const bool is_overlap =
          gjk_solver.shapeIntersect(tet1, tf1, tet2, tf2, &contacts);
      if (is_overlap) {
        collision_result[std::make_pair(i, j)] = std::move(contacts);
      }
    }
  }
  return !collision_result.empty();
}

template <typename S>
int isOverlapBVH(const std::shared_ptr<const fcl::CollisionGeometry<S>>& bvh,
                 const Transform3<S>& tf1, const Transform3<S>& tf2) {
  fcl::CollisionObject<S> object_1(
      bvh, tf1, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionObject<S> object_2(
      bvh, tf2, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionRequest<S> request;
  request.setMaxContactCount(std::numeric_limits<std::size_t>::max());
  fcl::CollisionResult<S> result;
  fcl::collide(&object_1, &object_2, request, result);
  return result.numContacts();
}

template <typename S>
bool isOverlapBVH(const std::shared_ptr<const fcl::CollisionGeometry<S>>& bvh,
                  const Transform3<S>& tf1, const Transform3<S>& tf2,
                  fcl::CollisionResult<S>& result) {
  fcl::CollisionObject<S> object_1(
      bvh, tf1, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionObject<S> object_2(
      bvh, tf2, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionRequest<S> request(std::numeric_limits<std::size_t>::max());
  request.useDefaultPenetration();
  fcl::collide(&object_1, &object_2, request, result);
  return result.isCollision();
}

template <typename S>
bool isOverlapBVH(const std::shared_ptr<const fcl::CollisionGeometry<S>>& bvh1,
                  const std::shared_ptr<const fcl::CollisionGeometry<S>>& bvh2,
                  const Transform3<S>& tf1, const Transform3<S>& tf2,
                  fcl::CollisionResult<S>& result) {
  fcl::CollisionObject<S> object_1(
      bvh1, tf1, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionObject<S> object_2(
      bvh2, tf2, typename fcl::CollisionObject<S>::GeometryLocalAABBComputed());
  fcl::CollisionRequest<S> request(std::numeric_limits<std::size_t>::max());
  request.useDefaultPenetration();
  fcl::collide(&object_1, &object_2, request, result);
  return result.isCollision();
}

template <typename BV>
void testBoxBVHModelOfTetrahedronCollision() {
  using S = typename BV::S;
  const Box<S> box(0.2, 0.4, 0.8);
  auto box_mesh = MakeBoxBVHTetrahedronModel<BV>(box, 0.1);
  EXPECT_EQ(box_mesh.build_state(), BVH_BUILD_STATE_PROCESSED);

  std::vector<Tetrahedron<S>> tetrahedrons;
  tetrahedrons.reserve(box_mesh.num_simplex());

  // Add elements
  for (auto idx = 0; idx < box_mesh.num_simplex(); idx++) {
    const auto simplex = box_mesh.getSimplex(idx);
    tetrahedrons.emplace_back(simplex[0], simplex[1], simplex[2], simplex[3]);
  }
  const auto bvhPtr = std::make_shared<const BVHModel<BV>>(std::move(box_mesh));

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  const int test_n = 1e2;
  const S eps = 1e-6;

  for (int i = 0; i < test_n; i++) {
    Transform3<S> tf1;
    fcl::test::generateRandomTransform(xyz_extent, tf1);
    Transform3<S> tf2;
    fcl::test::generateRandomTransform(xyz_extent, tf2);
    const int collision_count1 = isOverlapBVH<S>(bvhPtr, tf1, tf2);
    const bool is_overlap_1 = collision_count1 > 0;
    const int collision_count2 = isOverlap(tetrahedrons, tf1, tf2);
    EXPECT_EQ(collision_count1, collision_count2);

    fcl::CollisionResult<S> result;
    const bool is_overlap_3 = isOverlapBVH<S>(bvhPtr, tf1, tf2, result);
    EXPECT_EQ(is_overlap_1, is_overlap_3);
    CollisionResultMap<S> collision_result_map;
    const bool is_overlap_4 =
        isOverlap(tetrahedrons, tf1, tf2, collision_result_map);
    EXPECT_EQ(is_overlap_1, is_overlap_4);
    EXPECT_EQ(result.numContacts(), collision_result_map.size());

    for (const Contact<S>& bvh_contact : result.getContacts()) {
      const auto key = std::make_pair<int, int>(bvh_contact.b1, bvh_contact.b2);
      const auto iter = collision_result_map.find(key);
      const bool found = iter != collision_result_map.end();
      EXPECT_TRUE(found);
      if (found) {
        CollisionPenetrationContactData<S> contact = iter->second;
        EXPECT_NEAR(contact[0].penetration_depth, bvh_contact.penetration_depth,
                    eps);
        if (bvh_contact.penetration_depth > 1e-4) {
          EXPECT_TRUE(
              CompareMatrices(contact[0].normal, bvh_contact.normal, eps));
          EXPECT_TRUE(CompareMatrices(contact[0].pos, bvh_contact.pos, eps));
        }
      }
    }
  }
}

GTEST_TEST(MakeBoxBvhMeshTest, testBoxBVHModelOfTetrahedronCollision) {
  testBoxBVHModelOfTetrahedronCollision<fcl::AABB<double>>();
  testBoxBVHModelOfTetrahedronCollision<fcl::OBB<double>>();
  testBoxBVHModelOfTetrahedronCollision<fcl::RSS<double>>();
  testBoxBVHModelOfTetrahedronCollision<fcl::kIOS<double>>();
  testBoxBVHModelOfTetrahedronCollision<fcl::OBBRSS<double>>();
}

template <typename BV>
void testCylinderBVHModelOfTetrahedronCollision() {
  using S = typename BV::S;
  const S height = 2;
  const S radius = 1;
  const S resolution_hint = 2.0;  // Hint to coarsest mesh.
  auto cylinder_mesh = MakeCylinderBVHTetrahedronModel<BV>(
      Cylinder<S>(radius, height), resolution_hint);
  EXPECT_EQ(cylinder_mesh.build_state(), BVH_BUILD_STATE_PROCESSED);

  std::vector<Tetrahedron<S>> tetrahedrons;
  tetrahedrons.reserve(cylinder_mesh.num_simplex());

  // Add elements
  for (auto idx = 0; idx < cylinder_mesh.num_simplex(); idx++) {
    const auto simplex = cylinder_mesh.getSimplex(idx);
    tetrahedrons.emplace_back(simplex[0], simplex[1], simplex[2], simplex[3]);
  }
  const auto bvhPtr =
      std::make_shared<const BVHModel<BV>>(std::move(cylinder_mesh));

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  const int test_n = 1e2;
  const S eps = 1e-6;

  for (int i = 0; i < test_n; i++) {
    Transform3<S> tf1;
    fcl::test::generateRandomTransform(xyz_extent, tf1);
    Transform3<S> tf2;
    fcl::test::generateRandomTransform(xyz_extent, tf2);
    const int collision_count1 = isOverlapBVH<S>(bvhPtr, tf1, tf2);
    const bool is_overlap_1 = collision_count1 > 0;

    const int collision_count2 = isOverlap(tetrahedrons, tf1, tf2);
    EXPECT_EQ(collision_count1, collision_count2);

    fcl::CollisionResult<S> result;
    const bool is_overlap_3 = isOverlapBVH<S>(bvhPtr, tf1, tf2, result);
    EXPECT_EQ(is_overlap_1, is_overlap_3);
    CollisionResultMap<S> collision_result_map;
    const bool is_overlap_4 =
        isOverlap(tetrahedrons, tf1, tf2, collision_result_map);
    EXPECT_EQ(is_overlap_1, is_overlap_4);
    EXPECT_EQ(result.numContacts(), collision_result_map.size());

    for (const Contact<S>& bvh_contact : result.getContacts()) {
      const auto key = std::make_pair<int, int>(bvh_contact.b1, bvh_contact.b2);
      const auto iter = collision_result_map.find(key);
      const bool found = iter != collision_result_map.end();
      EXPECT_TRUE(found);
      if (found) {
        CollisionPenetrationContactData<S> contact = iter->second;
        EXPECT_NEAR(contact[0].penetration_depth, bvh_contact.penetration_depth,
                    eps);
        if (bvh_contact.penetration_depth > 1e-4) {
          EXPECT_TRUE(
              CompareMatrices(contact[0].normal, bvh_contact.normal, eps));
          EXPECT_TRUE(CompareMatrices(contact[0].pos, bvh_contact.pos, eps));
        }
      }
    }
  }
}

GTEST_TEST(MakeCylinderBvhMeshTest,
           testCylinderBVHModelOfTetrahedronCollision) {
  testCylinderBVHModelOfTetrahedronCollision<fcl::OBB<double>>();
  testCylinderBVHModelOfTetrahedronCollision<fcl::RSS<double>>();
  testCylinderBVHModelOfTetrahedronCollision<fcl::kIOS<double>>();
  testCylinderBVHModelOfTetrahedronCollision<fcl::OBBRSS<double>>();
}

template <typename BV>
void testBoxCylinderBVHModelOfTetrahedronCollision() {
  using S = typename BV::S;

  const Box<S> box(0.2, 0.4, 0.8);
  auto box_mesh = MakeBoxBVHTetrahedronModel<BV>(box, 0.1);
  EXPECT_EQ(box_mesh.build_state(), BVH_BUILD_STATE_PROCESSED);

  const S height = 1;
  const S radius = 1;
  Cylinder<S> cylinder(radius, height);
  const S resolution_hint = 0.5;  // Hint to coarsest mesh.
  auto cylinder_mesh =
      MakeCylinderBVHTetrahedronModel<BV>(cylinder, resolution_hint);
  EXPECT_EQ(cylinder_mesh.build_state(), BVH_BUILD_STATE_PROCESSED);

  std::vector<Tetrahedron<S>> box_tetrahedrons;
  box_tetrahedrons.reserve(box_mesh.num_simplex());
  std::vector<Tetrahedron<S>> cylinder_tetrahedrons;
  cylinder_tetrahedrons.reserve(cylinder_mesh.num_simplex());

  // Add elements
  for (auto idx = 0; idx < box_mesh.num_simplex(); idx++) {
    const auto simplex = box_mesh.getSimplex(idx);
    box_tetrahedrons.emplace_back(simplex[0], simplex[1], simplex[2],
                                  simplex[3]);
  }
  for (auto idx = 0; idx < cylinder_mesh.num_simplex(); idx++) {
    const auto simplex = cylinder_mesh.getSimplex(idx);
    cylinder_tetrahedrons.emplace_back(simplex[0], simplex[1], simplex[2],
                                       simplex[3]);
  }
  const auto box_geometry = std::make_shared<const Box<S>>(box);
  const auto cylinder_geometry = std::make_shared<const Cylinder<S>>(cylinder);
  const auto box_bvh =
      std::make_shared<const BVHModel<BV>>(std::move(box_mesh));
  const auto cylinder_bvh =
      std::make_shared<const BVHModel<BV>>(std::move(cylinder_mesh));

  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  const int test_n = 1e2;
  const S eps = 1e-4;
  detail::GJKSolver<S> gjk_solver;
  detail::ShapePairIntersectSolver<S> shapePairIntersectSolver(&gjk_solver);
  for (int i = 0; i < test_n; i++) {
    Transform3<S> tf1;
    fcl::test::generateRandomTransform(xyz_extent, tf1);
    Transform3<S> tf2;
    fcl::test::generateRandomTransform(xyz_extent, tf2);

    fcl::CollisionResult<S> result;
    const bool is_overlap_1 =
        isOverlapBVH<S>(box_bvh, cylinder_bvh, tf1, tf2, result);
    CollisionResultMap<S> collision_result_map;
    const bool is_overlap_2 = isOverlap(box_tetrahedrons, cylinder_tetrahedrons,
                                        tf1, tf2, collision_result_map);
    EXPECT_EQ(is_overlap_1, is_overlap_2);
    EXPECT_EQ(result.numContacts(), collision_result_map.size());
    for (const Contact<S>& bvh_contact : result.getContacts()) {
      const auto key = std::make_pair<int, int>(bvh_contact.b1, bvh_contact.b2);
      const auto iter = collision_result_map.find(key);
      const bool found = iter != collision_result_map.end();
      EXPECT_TRUE(found);
      if (found) {
        CollisionPenetrationContactData<S> contact = iter->second;
        EXPECT_NEAR(contact[0].penetration_depth, bvh_contact.penetration_depth,
                    eps);
        if (contact[0].penetration_depth > 1e-4 &&
            bvh_contact.penetration_depth > 1e-4) {
          EXPECT_TRUE(
              CompareMatrices(contact[0].normal, bvh_contact.normal, eps));
          EXPECT_TRUE(CompareMatrices(contact[0].pos, bvh_contact.pos, eps));
        }
      }
    }

    CollisionRequest<S> request;
    CollisionResult<S> result3;
    shapePairIntersectSolver.ShapeIntersect(&box, tf1, &cylinder, tf2, request,
                                            result3);
    const bool is_overlap_3 = result3.isCollision();

    EXPECT_EQ(is_overlap_1, is_overlap_3);

    fcl::CollisionResult<S> result4;
    const bool is_overlap_4 =
        isOverlapBVH<S>(box_geometry, cylinder_bvh, tf1, tf2, result4);
    EXPECT_EQ(is_overlap_1, is_overlap_4);

    fcl::CollisionResult<S> result5;
    const bool is_overlap_5 =
        isOverlapBVH<S>(box_bvh, cylinder_geometry, tf1, tf2, result5);
    EXPECT_EQ(is_overlap_1, is_overlap_5);
  }
}

GTEST_TEST(MakeCylinderBvhMeshTest,
           testBoxCylinderBVHModelOfTetrahedronCollision) {
  testBoxCylinderBVHModelOfTetrahedronCollision<fcl::OBB<double>>();
  testBoxCylinderBVHModelOfTetrahedronCollision<fcl::kIOS<double>>();
  testBoxCylinderBVHModelOfTetrahedronCollision<fcl::OBBRSS<double>>();
}

}  // namespace test
}  // namespace fcl

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
