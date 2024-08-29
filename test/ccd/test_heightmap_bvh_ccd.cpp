//
// Created by Wei Gao on 2024/6/18.
//
#include <gtest/gtest.h>

#include <set>

#include "fcl/narrowphase/continuous_collision.h"
#include "fcl/narrowphase/detail/ccd/heightmap_ccd_solver.h"
#include "fcl/narrowphase/detail/ccd/shape_pair_ccd.h"
#include "fcl/narrowphase/detail/traversal/heightmap/heightmap_solver.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

namespace fcl {

std::shared_ptr<heightmap::PointCloud> generateRandomPointCloud(
    std::size_t numPoints) {
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
  test::loadOBJFile(TEST_RESOURCES_DIR "/env.obj", p1, t1);

  // Update on vertex to scale it to 1[m]
  const S scale_down = S(1.0 / 4000.0);
  for (std::size_t i = 0; i < p1.size(); i++) {
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
int bvhHeightMapTranslationalCollisionCompareWithPrimitive(
    std::shared_ptr<const BVHModel<BV>> bvh,
    std::shared_ptr<const heightmap::LayeredHeightMap<typename BV::S>>
        height_map,
    const Transform3<typename BV::S>& tf_bvh,
    const Transform3<typename BV::S>& tf_hm,
    const TranslationalDisplacement<typename BV::S>& hm_displacement) {
  using S = typename BV::S;
  detail::TranslationalDisplacementHeightMapSolver<S> solver;
  fcl::HeightMapCollisionGeometry<S> hm_geometry(height_map);
  hm_geometry.computeLocalAABB();

  // First checking
  ContinuousCollisionRequest<S> request;
  request.num_max_contacts = 10000000;

  ContinuousCollisionResult<S> result0;
  solver.RunHeightMapObbBVH(&hm_geometry, tf_hm, hm_displacement, bvh.get(),
                            tf_bvh, request, result0);

  // Reverse order with fcl
  TranslationalDisplacement<S> bvh_equiv_disp;
  bvh_equiv_disp.scalar_displacement = hm_displacement.scalar_displacement;
  bvh_equiv_disp.unit_axis_in_shape1 =
      tf_bvh.linear().transpose() *
      (tf_hm.linear() * (-hm_displacement.unit_axis_in_shape1));

  {
    ContinuousCollisionResult<S> result1;
    fcl::translational_ccd(bvh.get(), tf_bvh, bvh_equiv_disp, &hm_geometry,
                           tf_hm, request, result1);
    EXPECT_EQ(result0.num_contacts(), result1.num_contacts());
  }

  // Verity collision with provided primitives
  for (std::size_t i = 0; i < result0.num_contacts(); i++) {
    const auto& contact_i = result0.raw_contacts()[i];
    // Get pixel and box
    heightmap::Pixel pixel_i = heightmap::decodePixel(contact_i.b1);
    Box<S> pixel_i_box;
    Transform3<S> pixel_i_box_tf_in_hm;
    auto ok = height_map->bottom().pixelToBox(pixel_i, pixel_i_box,
                                              pixel_i_box_tf_in_hm);
    EXPECT_TRUE(ok);
    Transform3<S> pixel_i_box_tf = tf_hm * pixel_i_box_tf_in_hm;
    pixel_i_box.computeLocalAABB();

    // Get the primitive
    auto bvh_primitive_id = contact_i.b2;
    const MeshSimplex& tri_id = bvh->simplex_indices()[bvh_primitive_id];
    const Vector3<S>& p1 = bvh->getVertex(tri_id[0]);
    const Vector3<S>& p2 = bvh->getVertex(tri_id[1]);
    const Vector3<S>& p3 = bvh->getVertex(tri_id[2]);
    Simplex<S> simplex(p1, p2, p3);

    ContinuousCollisionResult<S> box_result;
    detail::ContinuousContactMeta<S> box_contact_meta;
    using ShapePairSolver = detail::ShapePairTranslationalCollisionSolver<S>;
    ShapePairSolver::RunShapeSimplex(pixel_i_box, pixel_i_box_tf,
                                     hm_displacement, simplex, tf_bvh,
                                     box_contact_meta, request, box_result);
    bool intersect = box_result.num_contacts() > 0;
    EXPECT_TRUE(intersect);
  }

  // Do naive checking
  std::size_t naive_count = 0;
  for (int primitive_id = 0; primitive_id < bvh->num_simplex();
       primitive_id++) {
    const MeshSimplex& tri_id = bvh->simplex_indices()[primitive_id];
    const Vector3<S>& p1 = bvh->getVertex(tri_id[0]);
    const Vector3<S>& p2 = bvh->getVertex(tri_id[1]);
    const Vector3<S>& p3 = bvh->getVertex(tri_id[2]);
    fcl::TriangleP<S> triangle_p(p1, p2, p3);
    triangle_p.computeLocalAABB();
    ContinuousCollisionResult<S> triangle_result;
    solver.RunShapeHeightMap(&triangle_p, tf_bvh, bvh_equiv_disp, &hm_geometry,
                             tf_hm, request, triangle_result);
    naive_count += triangle_result.num_contacts();
  }

  EXPECT_EQ(naive_count, result0.num_contacts());
  // std::cerr << "#of contacts naive " << naive_count << std::endl;
  // std::cerr << "#of contacts " << result0.numContacts() << std::endl;
  return result0.num_contacts();
}

template <typename BV>
void heightMapBoxBVH_TranslationCollisionCompareWithPrimitive(int test_n) {
  // The basic geometric
  using S = typename BV::S;
  auto box_bvh = generateBoxBVHModel<BV>();
  auto obj_bvh = generateBVH_FromObj<BV>();
  auto point_cloud_0 = generateRandomPointCloud(10000);

  // Test candidate
  using BVH_HeightMapPair =
      std::pair<std::shared_ptr<const BVHModel<BV>>,
                std::shared_ptr<const heightmap::LayeredHeightMap<S>>>;
  std::vector<BVH_HeightMapPair> bvh_hm_pairs;

  // first one w box
  {
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(box_bvh, hm));
  }

  // second one w box
  {
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 0.004,
                                                               512, 256);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(box_bvh, hm));
  }

  // first one w obj
  {
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 512);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(obj_bvh, hm));
  }

  // second one w obj
  {
    auto hm = std::make_shared<heightmap::LayeredHeightMap<S>>(0.002, 0.004,
                                                               512, 256);
    hm->resetHeights();
    hm->updateHeightsByPointCloud3D(*point_cloud_0);
    bvh_hm_pairs.emplace_back(std::make_pair(obj_bvh, hm));
  }

  std::array<S, 6> extent{-1, -1, -0.1, 1, 1, 0.1};
  for (int test_idx = 0; test_idx < test_n; test_idx++) {
    Eigen::Transform<S, 3, Eigen::Isometry> tf1;
    test::generateRandomTransform(extent, tf1);
    Eigen::Transform<S, 3, Eigen::Isometry> tf2;
    test::generateRandomTransform(extent, tf2);

    // Random displacement
    TranslationalDisplacement<S> hm2_displacement;
    {
      Transform3<S> translation_pose;
      test::generateRandomTransform(extent, translation_pose);
      Matrix3<S> axis_in_cols = translation_pose.linear().matrix();

      // Setup the axis
      hm2_displacement.unit_axis_in_shape1 = axis_in_cols.col(0);
      hm2_displacement.scalar_displacement =
          std::abs(translation_pose.translation()[0]);
    }

    for (std::size_t i = 0; i < bvh_hm_pairs.size(); i++) {
      const auto& test_pair_i = bvh_hm_pairs[i];
      auto bvh_i = test_pair_i.first;
      auto height_map_i = test_pair_i.second;
      bvhHeightMapTranslationalCollisionCompareWithPrimitive<BV>(
          bvh_i, height_map_i, tf1, tf2, hm2_displacement);
    }
  }
}

}  // namespace fcl

GTEST_TEST(HeightMapBVH_Test, PrimitiveTest) {
  fcl::heightMapBoxBVH_TranslationCollisionCompareWithPrimitive<
      fcl::OBB<float>>(5);
  fcl::heightMapBoxBVH_TranslationCollisionCompareWithPrimitive<
      fcl::OBB<double>>(5);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}