//
// Created by mech-mind_gw on 3/22/2024.
//
#include <gtest/gtest.h>

#include "fcl/geometry/octree2/octree.h"
#include "fcl/geometry/octree2/octree_prune.h"
#include "fcl/geometry/octree2/octree_visit.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace octree2 {

template <typename S>
void octreeMetaInfoTest() {
  Vector3<S> resolution(0.4, 0.4, 0.4);

  {
    Octree<S> tree(resolution, 2);
    EXPECT_EQ(tree.n_layers(), 3u);
    EXPECT_EQ(tree.inner_nodes().size(), 1u);
    const auto& layer_meta = tree.layer_metas();
    EXPECT_EQ(layer_meta[0].depth, 0u);
    EXPECT_EQ(layer_meta[0].full_shape, 1u);
    EXPECT_EQ(layer_meta[0].half_shape, 0u);
    EXPECT_NEAR((layer_meta[0].resolution_xyz - 4 * resolution).norm(), 0.0,
                1e-7);

    EXPECT_EQ(layer_meta[1].depth, 1u);
    EXPECT_EQ(layer_meta[1].full_shape, 2u);
    EXPECT_EQ(layer_meta[1].half_shape, 1u);
    EXPECT_NEAR((layer_meta[1].resolution_xyz - 2 * resolution).norm(), 0.0,
                1e-7);

    // Check inv of all layers
    for (const auto& layer : layer_meta) {
      for (auto i = 0; i < 3; i++)
        EXPECT_NEAR(layer.resolution_xyz[i] * layer.inv_resolution_xyz[i], 1.0,
                    1e-7);
    }

    // Check child relation
    OctreeVoxel voxel;
    Vector3<S> point = Vector3<S>(-0.1, -0.1, -0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    auto child = tree.computeChildIndex(point, 0);
    auto child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 0u);
    EXPECT_EQ(child, child_voxel);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 7u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(0.1, -0.1, -0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 1u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 6u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(-0.1, 0.1, -0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 2u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 5u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(0.1, 0.1, -0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 3u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 4u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(-0.1, -0.1, 0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 4u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 3u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(0.1, -0.1, 0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 5u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 2u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(-0.1, 0.1, 0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 6u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 1u);
    EXPECT_EQ(child, child_voxel);

    point = Vector3<S>(0.1, 0.1, 0.1);
    EXPECT_TRUE(tree.computeVoxelCoordinate(point, voxel));
    child = tree.computeChildIndex(point, 0);
    child_voxel = tree.computeChildIndex(voxel, 0);
    EXPECT_EQ(child, 7u);
    child = tree.computeChildIndex(point, 1);
    child_voxel = tree.computeChildIndex(voxel, 1);
    EXPECT_EQ(child, 0u);
    EXPECT_EQ(child, child_voxel);
  }
}

template <typename S>
void collectLeafAABB(const Octree<S>& tree, std::vector<AABB<S>>& leaf_aabb) {
  leaf_aabb.clear();
  auto visit_fn = [&leaf_aabb](const AABB<S>& aabb, std::uint8_t,
                               bool is_leaf) -> bool {
    if (!is_leaf) return false;
    leaf_aabb.push_back(aabb);
    return false;
  };

  // Start visit
  visitOctree<S>(tree, visit_fn);
}

template <typename S>
void testContainmentNaive(const Octree<S>& tree,
                          const std::vector<Vector3<S>>& points_inserted,
                          const Vector3<S>* expected_leaf_aabb_size) {
  // Collect all the aabb
  std::vector<AABB<S>> leaf_aabb;
  collectLeafAABB(tree, leaf_aabb);

  // Naive query, can be very slot
  std::vector<bool> leaf_aabb_used;
  leaf_aabb_used.resize(leaf_aabb.size());
  std::fill(leaf_aabb_used.begin(), leaf_aabb_used.end(), false);

  for (const auto& point : points_inserted) {
    bool contained_in_aabb = false;
    for (std::size_t i = 0; i < leaf_aabb.size(); i++) {
      const auto& aabb = leaf_aabb[i];
      if (aabb.contain(point)) {
        leaf_aabb_used[i] = true;
        contained_in_aabb = true;
        break;
      }
    }

    EXPECT_TRUE(contained_in_aabb);
  }

  // Every aabb must be used
  for (std::size_t i = 0; i < leaf_aabb_used.size(); i++) {
    EXPECT_TRUE(leaf_aabb_used[i]);
  }

  // Check size
  if (expected_leaf_aabb_size == nullptr) return;
  const Vector3<S> aabb_expected_size = *expected_leaf_aabb_size;
  for (const auto& aabb : leaf_aabb) {
    const Vector3<S> this_aabb_size = aabb.max_ - aabb.min_;
    const auto size_diff = (this_aabb_size - aabb_expected_size).norm();
    EXPECT_NEAR(size_diff, 0, 1e-4);
  }
}

template <typename S>
void randomOccupiedTest(std::uint16_t bottom_half_shape = 2,
                        std::size_t test_n = 10) {
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  Vector3<S> resolution(scalar_resolution, scalar_resolution,
                        scalar_resolution);
  Octree<S> tree(resolution, bottom_half_shape);
  std::vector<Vector3<S>> points_inserted;
  for (std::size_t i = 0; i < test_n; i++) {
    Vector3<S> point_i;
    point_i.setRandom();
    point_i *= (0.99 * bottom_half_size);
    bool inserted_i = tree.Test_insertPointIntoTree(point_i);
    EXPECT_TRUE(inserted_i);
    EXPECT_TRUE(tree.isPointOccupied(point_i));
    points_inserted.push_back(point_i);
  }

  for (const auto& point_i : points_inserted) {
    EXPECT_TRUE(tree.isPointOccupied(point_i));
  }

  // Try copy
  {
    Octree<S> copied(tree);
    for (const auto& point_i : points_inserted) {
      EXPECT_TRUE(copied.isPointOccupied(point_i));
    }
  }

  const auto& raw_tree_nodes = tree.inner_nodes();
  const auto& raw_leaf_nodes = tree.leaf_nodes();
  std::cout << "Tree (inner, leaf) nodes size: (" << raw_tree_nodes.size()
            << ", " << raw_leaf_nodes.size() << ")" << std::endl;

  if (test_n < 1200) {
    testContainmentNaive<S>(tree, points_inserted, &resolution);
  }
}

template <typename S>
void randomRebuildTest(std::uint16_t bottom_half_shape = 2,
                       std::size_t test_n = 10) {
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  Vector3<S> resolution(scalar_resolution, scalar_resolution,
                        scalar_resolution);
  std::vector<Vector3<S>> points_inserted;
  for (std::size_t i = 0; i < test_n; i++) {
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
  Octree<S> tree(resolution, bottom_half_shape);
  tree.rebuildTree(point_fn, test_n);

  // They are all inserted
  for (const auto& point_i : points_inserted) {
    EXPECT_TRUE(tree.isPointOccupied(point_i));
  }

  if (test_n < 1200) {
    testContainmentNaive<S>(tree, points_inserted, nullptr);
  }
}

template <typename S>
void randomPruneTest(std::uint16_t bottom_half_shape = 2,
                     std::size_t test_n = 10) {
  const S scalar_resolution = 0.4;
  const S bottom_half_size = scalar_resolution * bottom_half_shape;
  Vector3<S> resolution(scalar_resolution, scalar_resolution,
                        scalar_resolution);
  std::vector<Vector3<S>> points_inserted;
  for (std::size_t i = 0; i < test_n; i++) {
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
  Octree<S> tree(resolution, bottom_half_shape);
  tree.rebuildTree(point_fn, test_n);

  // They are all inserted
  for (const auto& point_i : points_inserted) {
    EXPECT_TRUE(tree.isPointOccupied(point_i));
  }

  // Make the prune box
  AABB<S> obb_from_aabb;
  obb_from_aabb.min_.setConstant(-scalar_resolution * 0.3 * bottom_half_shape);
  obb_from_aabb.max_.setConstant(scalar_resolution * 0.3 * bottom_half_shape);
  Transform3<S> obb_pose;
  obb_pose.setIdentity();
  OBB<S> prune_obb;
  fcl::convertBV(obb_from_aabb, obb_pose, prune_obb);

  // Try prune
  OctreePruneInfo prune_info;
  pruneOctreeByOBB(tree, prune_obb, prune_info);

  auto visit_fn = [&prune_obb](const AABB<S>& aabb, std::uint8_t,
                               bool is_leaf) -> bool {
    if (!is_leaf) return false;
    const Vector3<S> center = aabb.center();
    EXPECT_FALSE(prune_obb.contain(center));
    return false;
  };
  visitOctree<S>(tree, &prune_info, visit_fn);

  // Try rebuild
  tree.rebuildAccordingToPruneInfo(prune_info);
  visitOctree<S>(tree, visit_fn);
}

}  // namespace octree2
}  // namespace fcl

GTEST_TEST(Octree2_ConstructByHandTest, MetaInfoTest) {
  fcl::octree2::octreeMetaInfoTest<float>();
  fcl::octree2::octreeMetaInfoTest<double>();
}

GTEST_TEST(Octree2_ConstructByHandTest, RandomInsertTest) {
  // Remove to avoid merge children in AABB size test
  // fcl::octree2::randomOccupiedTest<float>(2, 100);
  // fcl::octree2::randomOccupiedTest<double>(2, 100);
  fcl::octree2::randomOccupiedTest<float>(4, 100);
  fcl::octree2::randomOccupiedTest<double>(4, 100);
  fcl::octree2::randomOccupiedTest<float>(1024, 1000);
  fcl::octree2::randomOccupiedTest<double>(1024, 1000);
  fcl::octree2::randomOccupiedTest<float>(1024, 1000 * 100);
  fcl::octree2::randomOccupiedTest<double>(1024, 1000 * 100);
}

GTEST_TEST(Octree2_ConstructByHandTest, RandomRebuildTest) {
  fcl::octree2::randomRebuildTest<float>(2, 100);
  fcl::octree2::randomRebuildTest<double>(2, 100);
  fcl::octree2::randomRebuildTest<float>(4, 100);
  fcl::octree2::randomRebuildTest<double>(4, 100);
  fcl::octree2::randomRebuildTest<float>(1024, 1000);
  fcl::octree2::randomRebuildTest<double>(1024, 1000);
  fcl::octree2::randomRebuildTest<float>(1024, 1000 * 100);
  fcl::octree2::randomRebuildTest<double>(1024, 1000 * 100);
}

GTEST_TEST(Octree2_ConstructByHandTest, RandomPruneTest) {
  fcl::octree2::randomPruneTest<float>(2, 100);
  fcl::octree2::randomPruneTest<double>(2, 100);
  fcl::octree2::randomPruneTest<float>(4, 100);
  fcl::octree2::randomPruneTest<double>(4, 100);
  fcl::octree2::randomPruneTest<float>(1024, 1000);
  fcl::octree2::randomPruneTest<double>(1024, 1000);
  fcl::octree2::randomPruneTest<float>(1024, 1000 * 100);
  fcl::octree2::randomPruneTest<double>(1024, 1000 * 100);
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
