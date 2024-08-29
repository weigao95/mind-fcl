//
// Created by mech-mind_gw on 10/16/2023.
//
#include <gtest/gtest.h>
#include <unordered_map>
#include <iostream>

#include "eigen_matrix_compare.h"
#include "fcl/geometry/bvh/BVH_model.h"
#include "test_fcl_utility.h"

namespace fcl {

struct pair_hash {
  std::size_t operator()(const std::pair<int, int>& key) const {
    std::size_t a = std::hash<int>()(key.first);
    std::size_t b = std::hash<int>()(key.second);
    return a >= b ? a * a + a + b : a + b * b;
  }
};

template <typename BV>
void inspectBVHModelContainment(const BVHModel<BV>& bvh) {
  using S = typename BV::S;
  constexpr S containment_tol = 1e-4;
  const auto& primitive_indices = bvh.Test_primitiveIndices();

  auto test_simplex_containment = [](const BV& raw_bv,
                                     const Simplex<S>& simplex) -> void {
    for (auto j = 0; j < simplex.get_num_points(); j++) {
      const auto& point_j = simplex[j];
      const bool point_j_contained = raw_bv.contain(point_j);
      EXPECT_TRUE(point_j_contained);
    }
  };

  for (auto i = 0; i < bvh.getNumBVs(); i++) {
    const BVNode<BV>& bv_node_i = bvh.getBV(i);
    const auto& raw_bv_i = bv_node_i.bv;
    if (bv_node_i.isLeaf()) {
      const auto primitive_id = bv_node_i.primitiveId();
      const auto& raw_simplex = bvh.getSimplex(primitive_id);
      Simplex<S> simplex = raw_simplex;
      simplex.scaleWrtCenter(1.0 - containment_tol);
      test_simplex_containment(raw_bv_i, simplex);
    } else {
      const auto first_primitive = bv_node_i.first_primitive;
      const auto num_primitives = bv_node_i.num_primitives;
      for (auto k = 0; k < num_primitives; k++) {
        const auto raw_index_k = primitive_indices[first_primitive + k];
        const auto& raw_simplex = bvh.getSimplex(raw_index_k);
        Simplex<S> simplex = raw_simplex;
        simplex.scaleWrtCenter(1.0 - containment_tol);
        test_simplex_containment(raw_bv_i, simplex);
      }
    }
  }
}

template <typename BV>
void testBVHModelOneTetrahedron() {
  using S = typename BV::S;
  Vector3<S> a{0, 0, 0};
  Vector3<S> b{1, 0, 0};
  Vector3<S> c{0, 0.5, 0};
  Vector3<S> d{0, 0, 1.5};
  BVHModel<BV> bvh;

  // Start building
  int result = bvh.beginModel();
  EXPECT_EQ(result, BVH_OK);

  // Add one element
  result = bvh.addTetrahedron(a, b, c, d);
  EXPECT_EQ(result, BVH_OK);

  result = bvh.endModel();
  EXPECT_EQ(result, BVH_OK);

  // Should be one bv node
  EXPECT_EQ(bvh.num_simplex(), 1);
  EXPECT_EQ(bvh.getNumBVs(), 1);
  EXPECT_TRUE(bvh.getBV(0).isLeaf());
  inspectBVHModelContainment<BV>(bvh);
}

template <typename BV>
void testBVHModelTwoTetrahedron() {
  using S = typename BV::S;
  Vector3<S> a{0, 0, 0};
  Vector3<S> b{1, 0, 0};
  Vector3<S> c{0, 0.5, 0};
  Vector3<S> d{0, 0, 1.5};
  BVHModel<BV> bvh;

  // Start building
  int result = bvh.beginModel();
  EXPECT_EQ(result, BVH_OK);

  // Add one element
  result = bvh.addTetrahedron(a, b, c, d);
  EXPECT_EQ(result, BVH_OK);

  // Add another
  Vector3<S> displacement{2, 0, 0};
  result = bvh.addTetrahedron(a + displacement, b + displacement,
                              c + displacement, d + displacement);
  EXPECT_EQ(result, BVH_OK);

  result = bvh.endModel();
  EXPECT_EQ(result, BVH_OK);

  // Should be one bv node
  EXPECT_EQ(bvh.num_simplex(), 2);
  EXPECT_EQ(bvh.getNumBVs(), 3);
  inspectBVHModelContainment<BV>(bvh);

  // Check bvh topology
  const auto& bv_root = bvh.getBV(0);
  EXPECT_TRUE(!bv_root.isLeaf());
  EXPECT_EQ(bv_root.first_child, 1);
  EXPECT_EQ(bv_root.first_primitive, 0);
  EXPECT_EQ(bv_root.num_primitives, 2);

  EXPECT_TRUE(bvh.getBV(1).isLeaf());
  EXPECT_TRUE(bvh.getBV(2).isLeaf());
}

template <typename BV>
void testBVHModelArrayOfTetrahedron(int n_elem_x_axis = 10,
                                    int n_elem_y_axis = 10) {
  using S = typename BV::S;
  Vector3<S> a{0, 0, 0};
  Vector3<S> b{1, 0, 0};
  Vector3<S> c{0, 0.5, 0};
  Vector3<S> d{0, 0, 1.5};
  const S x_offset = 1.0;
  const S y_offset = 2.0;

  // Start building
  BVHModel<BV> bvh;
  int result = bvh.beginModel();
  EXPECT_EQ(result, BVH_OK);

  // Add elements
  for (auto x_idx = 0; x_idx < n_elem_x_axis; x_idx++) {
    for (auto y_idx = 0; y_idx < n_elem_y_axis; y_idx++) {
      Vector3<S> displacement{x_offset * S(x_idx), y_offset * S(y_idx), 0};
      result = bvh.addTetrahedron(a + displacement, b + displacement,
                                  c + displacement, d + displacement);
      EXPECT_EQ(result, BVH_OK);
    }
  }

  // Build the model
  result = bvh.endModel();
  EXPECT_EQ(result, BVH_OK);
  EXPECT_EQ(bvh.num_simplex(), n_elem_x_axis * n_elem_y_axis);
  EXPECT_EQ(bvh.getNumBVs(), n_elem_x_axis * n_elem_y_axis * 2 - 1);
  inspectBVHModelContainment<BV>(bvh);
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

template <typename S>
using CollisionResultMap =
    std::unordered_map<std::pair<int, int>, CollisionPenetrationContactData<S>,
                       pair_hash>;
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

template <typename BV>
void testBVHModelArrayOfTetrahedronCollision(int n_elem_x_axis = 10,
                                             int n_elem_y_axis = 10) {
  using S = typename BV::S;
  Vector3<S> a{0, 0, 0};
  Vector3<S> b{1, 0, 0};
  Vector3<S> c{0, 0.5, 0};
  Vector3<S> d{0, 0, 1.5};
  const S x_offset = 1.0;
  const S y_offset = 2.0;

  // Start building
  BVHModel<BV> bvh;
  std::vector<Tetrahedron<S>> tetrahedrons;
  tetrahedrons.reserve(n_elem_x_axis * n_elem_y_axis);
  int result = bvh.beginModel();
  EXPECT_EQ(result, BVH_OK);

  // Add elements
  for (auto x_idx = 0; x_idx < n_elem_x_axis; x_idx++) {
    for (auto y_idx = 0; y_idx < n_elem_y_axis; y_idx++) {
      Vector3<S> displacement{x_offset * S(x_idx), y_offset * S(y_idx), 0};
      result = bvh.addTetrahedron(a + displacement, b + displacement,
                                  c + displacement, d + displacement);
      tetrahedrons.emplace_back(a + displacement, b + displacement,
                                c + displacement, d + displacement);
      EXPECT_EQ(result, BVH_OK);
    }
  }

  // Build the model
  result = bvh.endModel();
  EXPECT_EQ(result, BVH_OK);

  const auto bvhPtr = std::make_shared<const BVHModel<BV>>(std::move(bvh));
  const std::array<S, 6> xyz_extent{-2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  const int test_n = 1e2;
  const S eps = 1e-6;

  for (int i = 0; i < test_n; i++) {
    Transform3<S> tf1;
    fcl::test::generateRandomTransform(xyz_extent, tf1);
    Transform3<S> tf2;
    fcl::test::generateRandomTransform(xyz_extent, tf2);
    const int collision_count1 = isOverlapBVH<S>(bvhPtr, tf1, tf2);
    const int collision_count2 = isOverlap(tetrahedrons, tf1, tf2);
    EXPECT_EQ(collision_count1, collision_count2);

    fcl::CollisionResult<S> result;
    const bool is_overlap_3 = isOverlapBVH<S>(bvhPtr, tf1, tf2, result);
    EXPECT_EQ(collision_count1, static_cast<int>(result.numContacts()));
    CollisionResultMap<S> collision_result_map;
    const bool is_overlap_4 =
        isOverlap(tetrahedrons, tf1, tf2, collision_result_map);
    EXPECT_EQ(result.numContacts(), collision_result_map.size());
    EXPECT_EQ(is_overlap_3, is_overlap_4);

    for (const Contact<S>& bvh_contact : result.getContacts()) {
      const auto key = std::make_pair<int, int>(bvh_contact.b1, bvh_contact.b2);
      const auto iter = collision_result_map.find(key);
      const bool find = iter != collision_result_map.end();
      EXPECT_TRUE(find);
      if (find) {
        CollisionPenetrationContactData<S> contact = iter->second;
        EXPECT_TRUE(
            CompareMatrices(contact[0].normal, bvh_contact.normal, eps));
        EXPECT_TRUE(CompareMatrices(contact[0].pos, bvh_contact.pos, eps));
        EXPECT_NEAR(contact[0].penetration_depth, bvh_contact.penetration_depth,
                    eps);
      }
    }
  }
}

}  // namespace fcl

GTEST_TEST(FCL_BVH_MODELS_TETRAHEDRON, building_one_tetrahedron) {
  fcl::testBVHModelOneTetrahedron<fcl::AABB<float>>();
  fcl::testBVHModelOneTetrahedron<fcl::OBB<float>>();
  fcl::testBVHModelOneTetrahedron<fcl::RSS<float>>();
  fcl::testBVHModelOneTetrahedron<fcl::kIOS<float>>();
  fcl::testBVHModelOneTetrahedron<fcl::OBBRSS<float>>();

  fcl::testBVHModelOneTetrahedron<fcl::AABB<double>>();
  fcl::testBVHModelOneTetrahedron<fcl::OBB<double>>();
  fcl::testBVHModelOneTetrahedron<fcl::RSS<double>>();
  fcl::testBVHModelOneTetrahedron<fcl::kIOS<double>>();
  fcl::testBVHModelOneTetrahedron<fcl::OBBRSS<double>>();
}

GTEST_TEST(FCL_BVH_MODELS_TETRAHEDRON, building_two_tetrahedron) {
  fcl::testBVHModelTwoTetrahedron<fcl::AABB<float>>();
  fcl::testBVHModelTwoTetrahedron<fcl::OBB<float>>();
  fcl::testBVHModelTwoTetrahedron<fcl::RSS<float>>();
  fcl::testBVHModelTwoTetrahedron<fcl::kIOS<float>>();
  fcl::testBVHModelTwoTetrahedron<fcl::OBBRSS<float>>();

  fcl::testBVHModelTwoTetrahedron<fcl::AABB<double>>();
  fcl::testBVHModelTwoTetrahedron<fcl::OBB<double>>();
  fcl::testBVHModelTwoTetrahedron<fcl::RSS<double>>();
  fcl::testBVHModelTwoTetrahedron<fcl::kIOS<double>>();
  fcl::testBVHModelTwoTetrahedron<fcl::OBBRSS<double>>();
}

GTEST_TEST(FCL_BVH_MODELS_TETRAHEDRON, building_array_tetrahedron) {
  fcl::testBVHModelArrayOfTetrahedron<fcl::AABB<double>>();
  fcl::testBVHModelArrayOfTetrahedron<fcl::OBB<double>>();
  fcl::testBVHModelArrayOfTetrahedron<fcl::RSS<double>>();
  fcl::testBVHModelArrayOfTetrahedron<fcl::kIOS<double>>();
  fcl::testBVHModelArrayOfTetrahedron<fcl::OBBRSS<double>>();
}

GTEST_TEST(FCL_BVH_MODELS_TETRAHEDRON,
           test_BVHModelArrayOfTetrahedron_Collision) {
  fcl::testBVHModelArrayOfTetrahedronCollision<fcl::AABB<double>>();
  fcl::testBVHModelArrayOfTetrahedronCollision<fcl::OBB<double>>();
  fcl::testBVHModelArrayOfTetrahedronCollision<fcl::RSS<double>>();
  fcl::testBVHModelArrayOfTetrahedronCollision<fcl::kIOS<double>>();
  fcl::testBVHModelArrayOfTetrahedronCollision<fcl::OBBRSS<double>>();
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
