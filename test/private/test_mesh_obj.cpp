//
// Created by mech-mind_gw on 7/26/2022.
//

#include <gtest/gtest.h>

#include <array>

#include "fcl/geometry/shape/convex.h"
#include "fcl/narrowphase/collision.h"
#include "test_fcl_utility.h"

template <typename T>
bool test_mesh_obj_collision() {
  // Get the path
  std::string test_model_dir;
  bool has_model = fcl::test::getTestModelDirectory(test_model_dir);
  if (!has_model) return -1;

  // Build bvh
  auto const bvhm = std::make_shared<fcl::BVHModel<fcl::OBB<T>>>();
  std::vector<fcl::Vector3<T>> mesh_points;
  std::vector<fcl::MeshSimplex> mesh_triangles;
  {
    // Get the model path
    std::string mesh_path = test_model_dir + "/4.STL";
    fcl::test::loadSTLMesh(mesh_path, mesh_points, mesh_triangles);
    std::cout << "# mesh points " << mesh_points.size() << std::endl;
    std::cout << "# mesh triangles " << mesh_triangles.size() << std::endl;

    bvhm->beginModel();
    bvhm->addSubModel(mesh_points, mesh_triangles);
    bvhm->endModel();
    bvhm->computeLocalAABB();
  }

  // Load objs
  std::shared_ptr<fcl::Convex<T>> cylinder_convex;
  {
    // Get the model path
    std::string mesh_path = test_model_dir + "/cylinder.obj";
    std::vector<fcl::Vector3f> points;
    std::vector<fcl::MeshSimplex> triangles;
    fcl::test::loadOBJFile(mesh_path.c_str(), points, triangles);
    std::cout << "# obj points " << points.size() << std::endl;
    std::cout << "# obj triangles " << triangles.size() << std::endl;

    // The format used by fcl::Convex
    std::shared_ptr<std::vector<fcl::Vector3<T>>> vertices =
        std::make_shared<std::vector<fcl::Vector3<T>>>();
    for (std::size_t i = 0; i < points.size(); i++) {
      vertices->push_back(
          fcl::Vector3<T>(points[i][0], points[i][1], points[i][2]));
    }

    int num_faces = triangles.size();
    std::shared_ptr<std::vector<int>> faces =
        std::make_shared<std::vector<int>>();
    for (std::size_t i = 0; i < triangles.size(); i++) {
      faces->push_back(3);
      faces->push_back(triangles[i][0]);
      faces->push_back(triangles[i][1]);
      faces->push_back(triangles[i][2]);
    }

    // Load it
    cylinder_convex =
        std::make_shared<fcl::Convex<T>>(vertices, num_faces, faces, true);
    cylinder_convex->computeLocalAABB();
  }

  // Original shape
  auto cylinder_original = std::make_shared<fcl::Cylinder<T>>(0.05f, 0.4f);
  cylinder_original->computeLocalAABB();

  // Transform for mesh
  fcl::Transform3<T> transform_mesh;
  transform_mesh.setIdentity();
  {
    transform_mesh.translation() =
        Eigen::Vector3<T>(0.56307351589202880859, -0.10914999991655349731,
                              0.19692927598953247070);
    Eigen::Quaternion<T> mesh_rotation(
        0.44619779485385335782, 0.00000000000000000000, 0.89493436767170619905,
        0.00000000000000000000);
    transform_mesh.linear() = mesh_rotation.toRotationMatrix();
  }

  // Transform for mesh
  fcl::Transform3<T> transform_cylinder;
  transform_cylinder.setIdentity();
  {
    transform_cylinder.translation() =
        Eigen::Vector3<T>(0.54610800000000003784, -0.14605199999999998739,
                              0.06137299999999999700);
  }

  // Run detection
  int n_contact_by_convex;
  {
    fcl::CollisionRequest<T> request;
    request.useDefaultPenetration();
    request.setMaxContactCount(INT32_MAX);
    fcl::CollisionResult<T> result;
    fcl::CollisionObject<T> mesh_object(bvhm, transform_mesh);
    fcl::CollisionObject<T> cylinder_object(cylinder_convex,
                                                transform_cylinder);
    fcl::collide(&mesh_object, &cylinder_object, request, result);
    n_contact_by_convex = result.numContacts();
  }

  // Do collision detection by primitive
  int n_contact_by_original = 0;
  {
    fcl::CollisionRequest<T> request;
    request.disablePenetration();
    request.setMaxContactCount(INT32_MAX);
    fcl::CollisionResult<T> result;
    fcl::CollisionObject<T> mesh_object(bvhm, transform_mesh);
    fcl::CollisionObject<T> cylinder_object(cylinder_original,
                                                transform_cylinder);
    fcl::collide(&mesh_object, &cylinder_object, request, result);
    n_contact_by_original = result.numContacts();
  }

  // Do collision detection by naive
  int n_contact_by_naive = 0;
  {
    fcl::detail::GJKSolver<T> solver;
    for (std::size_t i = 0; i < mesh_triangles.size(); i++) {
      const auto& triangle_i = mesh_triangles[i];
      auto intersect = solver.shapeTriangleIntersect(
          *cylinder_original.get(), transform_cylinder,
          mesh_points[triangle_i[0]], mesh_points[triangle_i[1]],
          mesh_points[triangle_i[2]], transform_mesh);
      if (intersect) {
        n_contact_by_naive += 1;
      }
    }
  }

  EXPECT_TRUE(n_contact_by_original == n_contact_by_naive);
  std::cout << "n_contact_by_original " << n_contact_by_original << std::endl;
  std::cout << "n_contact_by_convex " << n_contact_by_convex << std::endl;
  std::cout << "n_contact_by_naive " << n_contact_by_naive << std::endl;
  return 0;
}

GTEST_TEST(MeshConvexObjTest, FloatTest) {
  test_mesh_obj_collision<float>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
