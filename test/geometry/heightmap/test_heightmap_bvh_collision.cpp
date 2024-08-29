//
// Created by mech-mind_gw on 2/23/2023.
//

#include <gtest/gtest.h>
#include "fcl_resources/config.h"
#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "test_fcl_utility.h"
#include "fcl/narrowphase/collision.h"

namespace fcl {
namespace heightmap {

std::shared_ptr<PointCloud> generateRandomPointCloud(std::size_t numPoints) {
  auto cloud = std::make_shared<fcl::octomap::Pointcloud>();
  cloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    cloud->push_back(point[0], point[1], point[2] + 1);
  }
  return cloud;
}

template <typename BV>
std::shared_ptr<const BVHModel<BV>> generateBoxBVHModel() {
  using S = typename BV::S;
  Box<S> box(0.4, 0.4, 0.4);
  return test::generateBoxBVHModel<BV>(box);
}

template <typename BV>
std::shared_ptr<const BVHModel<BV>> generateBVH_FromObj() {
  using S = typename BV::S;
  std::vector<Vector3<S>> p1;
  std::vector<MeshSimplex> t1;
  test::loadOBJFile(TEST_RESOURCES_DIR"/env.obj", p1, t1);

  // Update on vertex to scale it to 1[m]
  const S scale_down = S(1.0/4000.0);
  for(std::size_t i = 0; i < p1.size(); i++) {
    p1[i] *= scale_down;
  }

  std::shared_ptr<BVHModel<BV>> model(new BVHModel<BV>());
  model->beginModel();
  model->addSubModel(p1, t1);
  model->endModel();
  model->computeLocalAABB();
  return model;
}

template <typename BV>
int bvhHeightMapCollisionCompareWithPrimitive(
    std::shared_ptr<const BVHModel<BV>> bvh,
    std::shared_ptr<const LayeredHeightMap<typename BV::S>> height_map,
    const Transform3<typename BV::S>& tf_bvh,
    const Transform3<typename BV::S>& tf_hm) {
  using S = typename BV::S;
  detail::GJKSolver<S> narrowphase_solver;
  detail::HeightMapCollisionSolver<S> solver(&narrowphase_solver);
  fcl::HeightMapCollisionGeometry<S> hm_geometry(height_map);
  hm_geometry.computeLocalAABB();

  // First checking
  CollisionRequest<S> request{UINT_MAX};
  CollisionResult<S> result0;
  solver.HeightMapBVHIntersect(&hm_geometry, bvh.get(), tf_hm, tf_bvh, request,
                               result0);

  // Reverse order with fcl
  CollisionResult<S> result1;
  fcl::collide(bvh.get(), tf_bvh, &hm_geometry, tf_hm, request, result1);
  EXPECT_EQ(result0.numContacts(), result1.numContacts());

  // Verity collision with provided primitives
  for (std::size_t i = 0; i < result0.numContacts(); i++) {
    const Contact<S>& contact_i = result0.getContact(i);
    // Get pixel and box
    Pixel pixel_i = decodePixel(contact_i.b1);
    Box<S> pixel_i_box;
    Transform3<S> pixel_i_box_tf_in_hm;
    auto ok = height_map->bottom().pixelToBox(pixel_i, pixel_i_box,
                                              pixel_i_box_tf_in_hm);
    EXPECT_TRUE(ok);
    Transform3<S> pixel_i_box_tf = tf_hm * pixel_i_box_tf_in_hm;

    // Get the primitive
    auto bvh_primitive_id = contact_i.b2;
    const MeshSimplex& tri_id = bvh->simplex_indices()[bvh_primitive_id];
    const Vector3<S>& p1 = bvh->getVertex(tri_id[0]);
    const Vector3<S>& p2 = bvh->getVertex(tri_id[1]);
    const Vector3<S>& p3 = bvh->getVertex(tri_id[2]);
    bool intersect = narrowphase_solver.shapeTriangleIntersect(
        pixel_i_box, pixel_i_box_tf, p1, p2, p3, tf_bvh);
    EXPECT_TRUE(intersect);
  }

  // Do naive checking
  std::size_t naive_count = 0;
  for (int primitive_id = 0; primitive_id < bvh->num_simplex(); primitive_id++) {
    const MeshSimplex& tri_id = bvh->simplex_indices()[primitive_id];
    const Vector3<S>& p1 = bvh->getVertex(tri_id[0]);
    const Vector3<S>& p2 = bvh->getVertex(tri_id[1]);
    const Vector3<S>& p3 = bvh->getVertex(tri_id[2]);
    fcl::TriangleP<S> triangle_p(p1, p2, p3);
    result1.clear();
    solver.ShapeHeightMapIntersect(triangle_p, &hm_geometry, tf_bvh, tf_hm,
                                   request, result1);
    naive_count += result1.numContacts();
  }

  EXPECT_EQ(naive_count, result0.numContacts());
  // std::cerr << "#of contacts naive " << naive_count << std::endl;
  // std::cerr << "#of contacts " << result0.numContacts() << std::endl;
  return result0.numContacts();
}

template <typename BV>
void heightMapBoxBVH_CompareWithPrimitive() {
  // The basic geometric
  using S = typename BV::S;
  auto box_bvh = generateBoxBVHModel<BV>();
  auto obj_bvh = generateBVH_FromObj<BV>();
  auto point_cloud_0 = generateRandomPointCloud(10000);

  // Test candidate
  using BVH_HeightMapPair =
      std::pair<std::shared_ptr<const BVHModel<BV>>,
                std::shared_ptr<const LayeredHeightMap<S>>>;
  std::vector<BVH_HeightMapPair> bvh_hm_pairs;

  // first one w box
  {
    auto hm = std::make_shared<LayeredHeightMap<S>>(0.002, 512);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(box_bvh, hm));
  }

  // second one w box
  {
    auto hm = std::make_shared<LayeredHeightMap<S>>(0.002, 0.004, 512, 256);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(box_bvh, hm));
  }

  // first one w obj
  {
    auto hm = std::make_shared<LayeredHeightMap<S>>(0.002, 512);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(obj_bvh, hm));
  }

  // second one w obj
  {
    auto hm = std::make_shared<LayeredHeightMap<S>>(0.002, 0.004, 512, 256);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(obj_bvh, hm));
  }

  // Loop over pairs
  int n_loop_per_pair = 5;
  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (int test_idx = 0; test_idx < n_loop_per_pair; test_idx++) {
    Eigen::Transform<S, 3, Eigen::Isometry> tf1;
    test::generateRandomTransform(extent, tf1);
    Eigen::Transform<S, 3, Eigen::Isometry> tf2;
    test::generateRandomTransform(extent, tf2);

    for (std::size_t i = 0; i < bvh_hm_pairs.size(); i++) {
      const auto& test_pair_i = bvh_hm_pairs[i];
      auto bvh_i = test_pair_i.first;
      auto height_map_i = test_pair_i.second;
      bvhHeightMapCollisionCompareWithPrimitive<BV>(
          bvh_i, height_map_i, tf1, tf2);
    }
  }
}

}  // namespace heightmap
}  // namespace fcl

GTEST_TEST(HeightMapBVH_Test, PrimitiveTest) {
  fcl::heightmap::heightMapBoxBVH_CompareWithPrimitive<fcl::AABB<float>>();
  fcl::heightmap::heightMapBoxBVH_CompareWithPrimitive<fcl::AABB<double>>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}