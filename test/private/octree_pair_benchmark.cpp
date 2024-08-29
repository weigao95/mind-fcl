//
// Created by wei on 2/2/22.
//
#include <array>
#include <chrono>

#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

std::shared_ptr<fcl::Octree2CollisionGeometry<float>> buildOctree2(
    float voxel_resolution = 0.001, std::uint16_t bottom_half_shape = 1024) {
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return nullptr;

  // Forward to builder
  std::string pcl_xyz_path = test_model_dir + "/pcl.xyz";
  return fcl::test::buildOctree2(pcl_xyz_path, voxel_resolution,
                                 bottom_half_shape);
}

std::shared_ptr<fcl::Octree2CollisionGeometry<float>> buildOctree2Binvox(
    std::uint16_t bottom_half_shape = 1024) {
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return nullptr;

  // Forward to builder
  std::string pcl_xyz_path = test_model_dir + "/X13.binvox";
  auto raw_tree =
      fcl::test::loadBinvoxAsOctree2(pcl_xyz_path, 1.0, bottom_half_shape);
  auto tree = std::make_shared<fcl::Octree2CollisionGeometry<float>>(raw_tree);
  return tree;
}

template <typename S>
void octree2PairBenchmark_0() {
  S resolution = 0.001;
  S delta_z = 0.0109;
  std::uint16_t bottom_half_shape = 1024;
  int max_n_contacts = 500000;
  auto fcl_tree_0 = buildOctree2(resolution, bottom_half_shape);
  if (fcl_tree_0 == nullptr) return;
  auto fcl_tree_1 = buildOctree2Binvox(128);
  if (fcl_tree_1 == nullptr) return;
  {
    // Should be ok as it is const ptr inside?
    fcl_tree_0->computeLocalAABB();
    fcl_tree_1->computeLocalAABB();

    // Make the object
    Eigen::Transform<S, 3, Eigen::Isometry> pose_0;
    Eigen::Transform<S, 3, Eigen::Isometry> pose_1;
    pose_0.setIdentity();
    pose_1.setIdentity();
    pose_0.translation().z() += delta_z;  // 2 [cm]
    fcl::CollisionObject<S> tree_obj_0(fcl_tree_0, pose_0);
    fcl::CollisionObject<S> tree_obj_1(fcl_tree_1, pose_1);

    // Try collide
    fcl::CollisionRequest<S> request;
    request.setMaxContactCount(max_n_contacts);
    // request.useDefaultPenetration();
    fcl::CollisionResult<S> result;

    // Run detection
    int test_n = 40;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto test_idx = 0; test_idx < test_n; test_idx++) {
      result.clear();
      fcl::collide(&tree_obj_0, &tree_obj_1, request, result);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms_time =
        std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
            .count();
    std::cout << "Octree2 time is ms " << ms_time << " for " << test_n
              << " runs " << std::endl;
    std::cerr << "The #of contacts is " << result.numContacts() << std::endl;
  }
}

template <typename S>
void obbIntersectBenchmark(std::size_t test_n) {
  fcl::AABB<S> aabb1;
  aabb1.min_.setConstant(-0.1);
  aabb1.max_.setConstant(0.2);
  fcl::AABB<S> aabb2 = aabb1;

  // Random pose
  const S extent_scalar = 0.5;
  std::array<S, 6> extent{-extent_scalar, -extent_scalar, -extent_scalar,
                          extent_scalar,  extent_scalar,  extent_scalar};
  std::vector<fcl::Transform3<S>> aabb1_poses;
  for (std::size_t i = 0; i < test_n; i++) {
    fcl::Transform3<S> aabb1_pose;
    fcl::test::generateRandomTransform(extent, aabb1_pose);
    aabb1_poses.push_back(aabb1_pose);
  }

  // Try benchmark
  fcl::OBB<S> obb1, obb2;
  fcl::convertBV(aabb2, fcl::Transform3<S>::Identity(), obb2);
  auto start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < aabb1_poses.size(); i++) {
    const auto& aabb1_pose = aabb1_poses[i];
    fcl::convertBV(aabb1, aabb1_pose, obb1);
    const bool intersect_i = obb1.overlap(obb2);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms_time =
      std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
          .count();
  std::cout << "OBB overlap time is ms " << ms_time << " for " << test_n
            << " runs " << std::endl;
}

int main() {
  // obbIntersectBenchmark<float>(10'000'000);
  // octreePairBenchmark_0<double>();
  octree2PairBenchmark_0<float>();
  // double cannot be tested due to building as float
  // octree2PairBenchmark_0<double>();
  return 0;
}