#include <gtest/gtest.h>

#include "fcl/geometry/heightmap/heightmap_collision_geometry.h"
#include "fcl/narrowphase/colliding_volume.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

namespace fcl {

std::shared_ptr<const octomap::Pointcloud> generateRandomPointCloud(
    std::size_t numPoints) {
  auto pointCloud = std::make_shared<octomap::Pointcloud>();
  pointCloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    pointCloud->push_back(point[0], point[1], point[2] + 1);
  }
  return pointCloud;
}

std::shared_ptr<const octomap::Pointcloud>
generateEvenlyDistributedPointCloud() {
  auto pointCloud = std::make_shared<octomap::Pointcloud>();
  for (float x = -1; x < 1; x += 0.001) {
    for (float y = -1; y < 1; y += 0.001) {
      pointCloud->push_back(x, y, 0);
    }
  }
  return pointCloud;
}

std::shared_ptr<const octomap::OcTree> generateLayerOcTree(double resolution) {
  auto tree = std::make_shared<octomap::OcTree>(resolution);
  for (float x = -1; x < 1; x += 0.001) {
    for (float y = -1; y < 1; y += 0.001) {
      tree->TEST_updateNode(octomap::point3d(x, y, 0));
    }
  }
  tree->updateInnerOccupancy();
  return tree;
}

template <typename S>
void simpleCollidingVolumeTest(
    const std::shared_ptr<const CollisionGeometry<S>>& o1,
    const std::shared_ptr<const CollisionGeometry<S>>& o2) {
  Transform3<S> tf_1;
  tf_1.setIdentity();
  Transform3<S> tf_2 = tf_1;
  tf_2.translation() = Vector3<S>(0, 0, 0.5);
  for (std::size_t request_count : {1, 10, 100, 1000}) {
    CollidingVolumeRequest<S> volume_request1(
        request_count, std::numeric_limits<double>::max());
    CollidingVolumeResult<S> volume_result1;
    approximateCollisionVolume(o1.get(), tf_1, o2.get(), tf_2, volume_request1,
                               volume_result1);

    CollisionRequest<S> default_request(request_count);
    CollisionResult<S> default_result;
    collide(o1.get(), tf_1, o2.get(), tf_2, default_request, default_result);
    EXPECT_EQ(volume_result1.numContacts(), default_result.numContacts());

    S min_volume = 1;
    CollidingVolumeRequest<S> volume_request2(
        request_count, request_count * min_volume, min_volume);
    CollidingVolumeResult<S> volume_result2;
    approximateCollisionVolume(o1.get(), tf_1, o2.get(), tf_2, volume_request2,
                               volume_result2);
    EXPECT_EQ(volume_result2.totalVolume(),
              default_result.numContacts() * min_volume);
  }
}

template <typename S>
void simpleCollidingVolumeTest() {
  const auto octree1 =
      std::shared_ptr<const octomap::OcTree>(test::generateOcTree(0.002));
  auto tree_geometry1 = std::make_shared<OcTree<S>>(octree1);
  tree_geometry1->computeLocalAABB();

  const auto octree2 =
      std::shared_ptr<const octomap::OcTree>(test::generateOcTree(0.002));
  auto tree_geometry2 = std::make_shared<OcTree<S>>(octree2);
  tree_geometry2->computeLocalAABB();
  simpleCollidingVolumeTest<S>(tree_geometry1, tree_geometry2);
  simpleCollidingVolumeTest<S>(tree_geometry1, tree_geometry1);

  const auto points1 = generateRandomPointCloud(100000);
  auto heightMap1 =
      std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
  heightMap1->updateHeightsByPointCloud3D(*points1);
  auto heightmap_geometry1 =
      std::make_shared<fcl::HeightMapCollisionGeometry<S>>(heightMap1);
  heightmap_geometry1->computeLocalAABB();

  const auto points2 = generateRandomPointCloud(100000);
  auto heightMap2 =
      std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
  heightMap2->updateHeightsByPointCloud3D(*points2);
  auto heightmap_geometry2 =
      std::make_shared<fcl::HeightMapCollisionGeometry<S>>(heightMap2);
  heightmap_geometry2->computeLocalAABB();
  simpleCollidingVolumeTest<S>(heightmap_geometry1, heightmap_geometry2);
  simpleCollidingVolumeTest<S>(heightmap_geometry1, heightmap_geometry1);
}

template <typename S>
void boxPlainHeightMapCollidingVolumeTest() {
  auto heightMap =
      std::make_shared<heightmap::LayeredHeightMap<S>>(0.001, 1024);
  const S height_of_heightMap = 1.;
  auto update_heightmap_visitor =
      [&height_of_heightMap](const heightmap::Pixel& pixel,
                             const heightmap::Point2D<S>& box_bottom_center,
                             uint16_t old_height_in_mm,
                             uint16_t& new_height_in_mm) -> bool {
    (void) (pixel);
    (void) (box_bottom_center);
    (void) (old_height_in_mm);
    new_height_in_mm = height_of_heightMap * 1000;
    return false;
  };
  heightMap->updateHeightsByBottomLayerUpdateFunctor(update_heightmap_visitor);
  auto heightmap_geometry =
      std::make_shared<fcl::HeightMapCollisionGeometry<S>>(heightMap);
  heightmap_geometry->computeLocalAABB();

  std::vector<Vector3<S>> sides{
      Vector3<S>(0.5, 0.5, 0.5),
      Vector3<S>(0.3, 0.7, 0.5),
  };

  for (const auto& side : sides) {
    auto box = std::make_shared<Box<S>>(side);
    box->computeLocalAABB();
    for (const bool use_rotate_z : {true, false}) {
      for (int i = 0; i < 10; i++) {
        auto x = test::rand_interval(-0.5, 0.5);
        auto y = test::rand_interval(-0.5, 0.5);
        auto z = test::rand_interval(-0.5, 0.5);

        Vector3<S> T(x, y, z);
        Transform3<S> tf_box;
        tf_box.setIdentity();
        if (use_rotate_z) {
          Matrix3<S> R;
          auto c = test::rand_interval((S)0, 2 * constants<S>::pi());
          test::eulerToMatrix<S>(0, 0, c, R);
          tf_box.linear() = R;
        }
        tf_box.translation() = T;

        CollidingVolumeRequest<S> volume_request(
            std::numeric_limits<std::size_t>::max(),
            std::numeric_limits<double>::max());
        CollidingVolumeResult<S> volume_result1;
        approximateCollisionVolume(heightmap_geometry.get(),
                                   Transform3<S>::Identity(), box.get(), tf_box,
                                   volume_request, volume_result1);
        CollidingVolumeResult<S> volume_result2;
        approximateCollisionVolume(heightmap_geometry.get(),
                                   Transform3<S>::Identity(), box.get(), tf_box,
                                   volume_request, volume_result2);

        const AABB<S> aabb_after_transform =
            translate(box->aabb_local, tf_box.translation());
        S intersecting_part_height =
            std::min(height_of_heightMap, aabb_after_transform.max_[2]) -
            std::max((S)0, aabb_after_transform.min_[2]);
        intersecting_part_height = std::max(intersecting_part_height, (S)0.);

        const double volume =
            box->side[0] * box->side[1] * intersecting_part_height;

        EXPECT_TRUE(volume_result1.totalVolume() ==
                    volume_result2.totalVolume());
        EXPECT_TRUE((volume_result1.totalVolume() * (1 + 1e-3)) >= volume);
        if (!use_rotate_z) {
          EXPECT_TRUE((volume_result1.totalVolume() * (0.9)) <= volume);
        }
      }
    }
  }
}

template <typename S>
void boxPlainOctreeCollidingVolumeTest() {
  const double resolution = 0.001;
  auto layer_octree = generateLayerOcTree(resolution);
  auto tree = std::make_shared<OcTree<S>>(layer_octree);
  tree->computeLocalAABB();

  std::vector<Vector3<S>> sides{
      Vector3<S>(0.5, 0.5, 0.5),
      Vector3<S>(0.3, 0.7, 0.5),
  };

  for (const auto& side : sides) {
    auto box = std::make_shared<Box<S>>(side);
    box->computeLocalAABB();
    for (const bool use_rotate_z : {true, false}) {
      for (int i = 0; i < 10; i++) {
        auto x = test::rand_interval(-0.5, 0.5);
        auto y = test::rand_interval(-0.5, 0.5);

        Vector3<S> T(x, y, 0);
        Transform3<S> tf_box;
        tf_box.setIdentity();
        if (use_rotate_z) {
          Matrix3<S> R;
          auto c = test::rand_interval((S)0, 2 * constants<S>::pi());
          test::eulerToMatrix<S>(0, 0, c, R);
          tf_box.linear() = R;
        }
        tf_box.translation() = T;

        CollidingVolumeRequest<S> volume_request(
            std::numeric_limits<std::size_t>::max(),
            std::numeric_limits<double>::max());
        CollidingVolumeResult<S> volume_result1;
        approximateCollisionVolume(tree.get(), Transform3<S>::Identity(),
                                   box.get(), tf_box, volume_request,
                                   volume_result1);
        CollidingVolumeResult<S> volume_result2;
        approximateCollisionVolume(tree.get(), Transform3<S>::Identity(),
                                   box.get(), tf_box, volume_request,
                                   volume_result2);

        const double volume =
            box->side[0] * box->side[1] * static_cast<S>(resolution);

        EXPECT_TRUE(volume_result1.totalVolume() ==
                    volume_result2.totalVolume());
        EXPECT_TRUE((volume_result1.totalVolume() * (1 + 1e-3)) >= volume);
        if (!use_rotate_z) {
          EXPECT_TRUE((volume_result1.totalVolume() * (0.9)) <= volume);
        }
      }
    }
  }
}

}  // namespace fcl

GTEST_TEST(CollisionVolumeTest, SimpleCollidingVolumeTest) {
  fcl::simpleCollidingVolumeTest<float>();
  fcl::simpleCollidingVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, BoxPlainHeightMapCollidingVolumeTest) {
  fcl::boxPlainHeightMapCollidingVolumeTest<float>();
  fcl::boxPlainHeightMapCollidingVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, BoxPlainOctreeCollidingVolumeTest) {
  fcl::boxPlainOctreeCollidingVolumeTest<float>();
  fcl::boxPlainOctreeCollidingVolumeTest<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
