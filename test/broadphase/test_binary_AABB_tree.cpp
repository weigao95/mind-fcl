//
// Created by Wei Gao on 2024/8/10.
//
#include <gtest/gtest.h>

#include "fcl/broadphase/binary_AABB_tree.h"

namespace fcl {
namespace detail {

template <typename S>
void simpleObjectArrayTest(std::uint32_t n_objects) {
  std::vector<BroadphaseObjectInfo<S>> objects;
  Vector3<S> center_begin(0, 0, 0);
  Vector3<S> center_delta(-0.4, 0, 0);
  Vector3<S> half_size(0.1, 0.1, 0.1);
  for (std::uint32_t i = 0; i < n_objects; i++) {
    Vector3<S> center_i = center_begin + S(i) * center_delta;
    BroadphaseObjectInfo<S> object_i;
    object_i.bv.min_ = center_i - half_size;
    object_i.bv.max_ = center_i + half_size;
    object_i.user_id = i;
    objects.emplace_back(object_i);
  }

  BinaryAABB_Tree<S, SimpleVectorObjectAllocator> tree;
  tree.Rebuild(objects.data(), objects.size());
  EXPECT_TRUE(tree.SanityCheck());

  // Note that this is an mutated interface, objects should be sorted by x
  for (std::uint32_t i = 1; i < n_objects; i++) {
    EXPECT_LT(objects[i - 1].bv.center()[0], objects[i].bv.center()[0]);
  }

  // Inspect the tree
  auto visit_tree = [&center_delta, &center_begin](
                        const AABB<S>& node_aabb, bool is_leaf,
                        std::uint64_t user_id_if_leaf, bool& current_node_done,
                        bool& overall_done) -> void {
    current_node_done = false;
    overall_done = false;
    if (!is_leaf) return;

    // Leaf node, check center
    const Vector3<S> center = node_aabb.center();
    const Vector3<S> center_from_id =
        center_begin + S(user_id_if_leaf) * center_delta;
    EXPECT_NEAR(center.x(), center_from_id.x(), 1e-5);
    EXPECT_NEAR(center.y(), center_from_id.y(), 1e-5);
    EXPECT_NEAR(center.z(), center_from_id.z(), 1e-5);
  };
  tree.VisitTree(visit_tree);

  // Try collision with single objects
  auto single_object_collision_fn = [](std::uint64_t leaf1, std::uint64_t leaf2,
                                       void*) -> bool {
    const auto id1 = leaf1;
    const auto id2 = leaf2;
    EXPECT_EQ(id1 * id1 + 1, id2);
    return false;
  };

  for (std::uint32_t i = 0; i < n_objects; i++) {
    Vector3<S> center_i = center_begin + S(i) * center_delta;
    BroadphaseObjectInfo<S> object_i;
    object_i.bv.min_ = center_i - half_size;
    object_i.bv.max_ = center_i + half_size;
    object_i.user_id = i * i + 1;
    tree.SingleObjectCollision(object_i, single_object_collision_fn, nullptr);
  }

  // Try collision with tree
  auto tree_collision_fn = [](std::uint64_t leaf1, std::uint64_t leaf2,
                              void*) -> bool {
    const auto id1 = leaf1;
    const auto id2 = leaf2;
    EXPECT_EQ(id1, id2);
    return false;
  };
  tree.SelfCollision(tree_collision_fn, nullptr);
}

template <typename S>
void mutationTest(std::uint32_t n_test) {
  const auto n_objects = n_test * 2;
  std::vector<BroadphaseObjectInfo<S>> objects;
  Vector3<S> center_begin(0, 0, 0);
  Vector3<S> center_delta(-0.4, 0, 0);
  Vector3<S> half_size(0.1, 0.1, 0.1);
  for (std::uint32_t i = 0; i < n_objects; i++) {
    Vector3<S> center_i = center_begin + S(i) * center_delta;
    BroadphaseObjectInfo<S> object_i;
    object_i.bv.min_ = center_i - half_size;
    object_i.bv.max_ = center_i + half_size;
    object_i.user_id = i;
    objects.emplace_back(object_i);
  }

  // Build the tree with half of the objects
  const auto n_half_objects = n_test;
  BinaryAABB_Tree<S, SimpleVectorObjectAllocator> tree;
  tree.Rebuild(objects.data(), n_half_objects);
  EXPECT_TRUE(tree.SanityCheck());

  // Inspect the tree
  auto visit_tree = [&center_delta, &center_begin](
                        const AABB<S>& node_aabb, bool is_leaf,
                        std::uint64_t user_id_if_leaf, bool& current_node_done,
                        bool& overall_done) -> void {
    current_node_done = false;
    overall_done = false;
    if (!is_leaf) return;

    // Leaf node, check center
    const Vector3<S> center = node_aabb.center();
    const Vector3<S> center_from_id =
        center_begin + S(user_id_if_leaf) * center_delta;
    EXPECT_NEAR(center.x(), center_from_id.x(), 1e-5);
    EXPECT_NEAR(center.y(), center_from_id.y(), 1e-5);
    EXPECT_NEAR(center.z(), center_from_id.z(), 1e-5);
  };
  tree.VisitTree(visit_tree);

  using UpdateState =
      typename BinaryAABB_Tree<S, SimpleVectorObjectAllocator>::TreeUpdateState;

  // Add objects
  {
    UpdateState state;
    tree.PrepareAddNewObjects(objects.data() + n_half_objects, n_half_objects,
                              state);
    tree.ApplyAddNewObjects(state);
    tree.VisitTree(visit_tree);
    EXPECT_EQ(tree.n_leaves(), n_objects);
    EXPECT_TRUE(tree.SanityCheck());
  }

  // Remove object and update structure
  {
    EXPECT_TRUE(tree.RemoveObject(0));
    UpdateState state;
    tree.PrepareUpdateStructure(state);
    tree.ApplyUpdateStructure(std::move(state));
    EXPECT_EQ(tree.n_leaves(), n_objects - 1);
    EXPECT_TRUE(tree.SanityCheck());
    tree.VisitTree(visit_tree);
  }
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(BinaryAABB_TreeTest, SimpleTest) {
  fcl::detail::simpleObjectArrayTest<float>(1);
  fcl::detail::simpleObjectArrayTest<float>(2);
  fcl::detail::simpleObjectArrayTest<float>(3);
  fcl::detail::simpleObjectArrayTest<float>(1000);
  fcl::detail::simpleObjectArrayTest<float>(10000);
}

GTEST_TEST(BinaryAABB_TreeTest, MutationTest) {
  fcl::detail::mutationTest<float>(1);
  fcl::detail::mutationTest<float>(2);
  fcl::detail::mutationTest<float>(3);
  fcl::detail::mutationTest<float>(1000);
  fcl::detail::mutationTest<float>(10000);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}