//
// Created by Wei Gao on 2024/8/15.
//
#include "fcl/broadphase/broadphase_AABB_tree.h"
#include "fcl/geometry/geometric_shape_to_BVH_model.h"
#include "fcl/narrowphase/collision.h"
#include "legacy_broadphase_bruteforce.h"
#include "test_fcl_utility.h"

namespace fcl {
namespace test {

template <typename S>
void generateEnvironments(std::vector<CollisionObject<S>>& env, S env_scale,
                          std::size_t n) {
  S extents[] = {-env_scale, env_scale,  -env_scale,
                 env_scale,  -env_scale, env_scale};
  aligned_vector<Transform3<S>> transforms;

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Box<S>* box = new Box<S>(5, 10, 20);
    box->computeLocalAABB();
    env.push_back(CollisionObject<S>(std::shared_ptr<CollisionGeometry<S>>(box),
                                     transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Sphere<S>* sphere = new Sphere<S>(30);
    sphere->computeLocalAABB();
    env.push_back(CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(sphere), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  for (std::size_t i = 0; i < n; ++i) {
    Cylinder<S>* cylinder = new Cylinder<S>(10, 40);
    cylinder->computeLocalAABB();
    env.push_back(CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(cylinder), transforms[i]));
  }
}

template <typename S>
void generateEnvironmentsMesh(std::vector<CollisionObject<S>>& env, S env_scale,
                              std::size_t n) {
  S extents[] = {-env_scale, env_scale,  -env_scale,
                 env_scale,  -env_scale, env_scale};
  aligned_vector<Transform3<S>> transforms;

  generateRandomTransforms(extents, transforms, n);
  Box<S> box(5, 10, 20);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, box, Transform3<S>::Identity());
    model->computeLocalAABB();
    env.push_back(CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  Sphere<S> sphere(30);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, sphere, Transform3<S>::Identity(), 16, 16);
    model->computeLocalAABB();
    env.push_back(CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }

  generateRandomTransforms(extents, transforms, n);
  Cylinder<S> cylinder(10, 40);
  for (std::size_t i = 0; i < n; ++i) {
    BVHModel<OBBRSS<S>>* model = new BVHModel<OBBRSS<S>>();
    generateBVHModel(*model, cylinder, Transform3<S>::Identity(), 16, 16);
    model->computeLocalAABB();
    env.push_back(CollisionObject<S>(
        std::shared_ptr<CollisionGeometry<S>>(model), transforms[i]));
  }
}

}  // namespace test
}  // namespace fcl

namespace fcl {

template <typename S>
void broadphaseCollisionTestBinaryAABB_Tree(S env_scale, std::size_t env_size,
                                            std::size_t query_size,
                                            std::size_t num_max_contacts = 1,
                                            bool exhaustive = false,
                                            bool use_mesh = false) {
  std::vector<CollisionObject<S>> env;
  if (use_mesh)
    test::generateEnvironmentsMesh(env, env_scale, env_size);
  else
    test::generateEnvironments(env, env_scale, env_size);

  std::vector<CollisionObject<S>> query;
  if (use_mesh)
    test::generateEnvironmentsMesh(query, env_scale, query_size);
  else
    test::generateEnvironments(query, env_scale, query_size);

  // update the environment
  S delta_angle_max = 10 / 360.0 * 2 * constants<S>::pi();
  S delta_trans_max = 0.01 * env_scale;
  for (std::size_t i = 0; i < env.size(); ++i) {
    S rand_angle_x = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_x = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;
    S rand_angle_y = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_y = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;
    S rand_angle_z = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_z = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;

    Matrix3<S> dR(AngleAxis<S>(rand_angle_x, Vector3<S>::UnitX()) *
                  AngleAxis<S>(rand_angle_y, Vector3<S>::UnitY()) *
                  AngleAxis<S>(rand_angle_z, Vector3<S>::UnitZ()));
    Vector3<S> dT(rand_trans_x, rand_trans_y, rand_trans_z);

    Matrix3<S> R = env[i].getRotation();
    Vector3<S> T = env[i].getTranslation();
    env[i].setTransform(dR * R, dR * T + dT);
    env[i].computeAABB();
  }

  for (std::size_t i = 0; i < query.size(); ++i) {
    S rand_angle_x = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_x = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;
    S rand_angle_y = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_y = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;
    S rand_angle_z = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_angle_max;
    S rand_trans_z = 2 * (rand() / (S)RAND_MAX - 0.5) * delta_trans_max;

    Matrix3<S> dR(AngleAxis<S>(rand_angle_x, Vector3<S>::UnitX()) *
                  AngleAxis<S>(rand_angle_y, Vector3<S>::UnitY()) *
                  AngleAxis<S>(rand_angle_z, Vector3<S>::UnitZ()));
    Vector3<S> dT(rand_trans_x, rand_trans_y, rand_trans_z);

    Matrix3<S> R = env[i].getRotation();
    Vector3<S> T = env[i].getTranslation();
    query[i].setTransform(dR * R, dR * T + dT);
    query[i].computeAABB();
  }

  // For manager objects
  std::vector<BroadphaseObjectInfo<S>> broadphase_objects;
  for (std::size_t i = 0; i < env.size(); i++) {
    BroadphaseObjectInfo<S> object_i;
    object_i.bv = env[i].getAABB();
    object_i.user_id = i;
    broadphase_objects.emplace_back(object_i);
  }

  std::vector<BroadphaseObjectInfo<S>> broadphase_objects_query;
  for (std::size_t i = 0; i < query.size(); i++) {
    BroadphaseObjectInfo<S> object_i;
    object_i.bv = query[i].getAABB();
    object_i.user_id = i;
    broadphase_objects_query.emplace_back(object_i);
  }

  // Compare with dynamic manager
  std::vector<CollisionObject<S>*> env_ptr;
  for (std::size_t i = 0; i < env.size(); i++) {
    env_ptr.push_back(&env[i]);
  }

  std::vector<CollisionObject<S>*> query_ptr;
  for (std::size_t i = 0; i < query.size(); i++) {
    query_ptr.push_back(&query[i]);
  }

  NaiveCollisionManager<S> manager, query_manager;
  manager.registerObjects(env_ptr);
  query_manager.registerObjects(query_ptr);

  DefaultCollisionData<S> data_manager;
  {
    auto& data = data_manager;
    if (exhaustive)
      data.request.setMaxContactCount(100000);
    else
      data.request.setMaxContactCount(num_max_contacts);

    auto callback = [](CollisionObject<S>* o1, CollisionObject<S>* o2,
                       void* data) {
      assert(data != nullptr);
      auto* collision_data = static_cast<DefaultCollisionData<S>*>(data);
      const CollisionRequest<S>& request = collision_data->request;
      CollisionResult<S>& result = collision_data->result;

      if (collision_data->done) return true;
      collide(o1, o2, request, result);

      if (result.isCollision() &&
          result.numContacts() >= request.maxNumContacts()) {
        collision_data->done = true;
      }

      return collision_data->done;
    };

    test::Timer timer;
    timer.start();
    manager.collide(&query_manager, &data, callback);
    timer.stop();
    // std::cout << "#of contacts " << data.result.numContacts() << std::endl;
    std::cout << "Elapsed time for collision manager in micro "
              << timer.getElapsedTimeInMicroSec() << std::endl;
  }

  // Build the tree
  BroadphaseAABB_Tree<S> tree, query_tree;
  tree.Rebuild(broadphase_objects.data(), broadphase_objects.size());
  query_tree.Rebuild(broadphase_objects_query.data(),
                     broadphase_objects_query.size());

  // Make collision fn
  auto collision_fn = [&env, &query](std::uint64_t leaf1, std::uint64_t leaf2,
                                     void* collision_fn_data) -> bool {
    auto* collision_data =
        static_cast<DefaultCollisionData<S>*>(collision_fn_data);
    if (collision_data->done) return true;

    // Run collision
    const CollisionRequest<S>& request = collision_data->request;
    CollisionResult<S>& result = collision_data->result;
    const auto& o1 = env[leaf1];
    const auto& o2 = query[leaf2];
    collide(&o1, &o2, request, result);

    // Update result
    if (result.isCollision() &&
        result.numContacts() >= request.maxNumContacts()) {
      collision_data->done = true;
    }

    // Termination
    return collision_data->done;
  };

  // Checking for AABB tree
  DefaultCollisionData<S> data_AABB_tree;
  {
    auto& data = data_AABB_tree;
    if (exhaustive)
      data.request.setMaxContactCount(100000);
    else
      data.request.setMaxContactCount(num_max_contacts);

    test::Timer timer;
    timer.start();
    tree.TreeCollision(query_tree, collision_fn, &data);
    timer.stop();
    // std::cout << "#of contacts " << data.result.numContacts() << std::endl;
    std::cout << "Elapsed time for binary AABB tree in micro "
              << timer.getElapsedTimeInMicroSec() << std::endl;
  }

  // Check result
  EXPECT_EQ(data_manager.result.numContacts(),
            data_AABB_tree.result.numContacts());

  // Checking with naive query
  DefaultCollisionData<S> data_AABB_tree_single_object;
  {
    // Assign option
    auto& data = data_AABB_tree_single_object;
    if (exhaustive)
      data.request.setMaxContactCount(100000);
    else
      data.request.setMaxContactCount(num_max_contacts);

    // Query size
    test::Timer timer;
    timer.start();
    for (std::uint32_t i = 0; i < query.size(); i++) {
      tree.SingleObjectCollision(broadphase_objects_query[i], collision_fn,
                                 &data);
    }
    timer.stop();
    // std::cout << "#of contacts in single object " <<
    // data.result.numContacts()
    //           << std::endl;
    std::cout << "Elapsed time for binary AABB tree single object in micro "
              << timer.getElapsedTimeInMicroSec() << std::endl;
  }

  // Check result
  EXPECT_EQ(data_AABB_tree_single_object.result.numContacts(),
            data_AABB_tree.result.numContacts());
}

}  // namespace fcl

/// check the update, only return collision or not
GTEST_TEST(CollisionBinaryAABB_Tree,
           test_core_bf_broad_phase_update_collision_binary) {
#ifdef NDEBUG
  fcl::broadphaseCollisionTestBinaryAABB_Tree<double>(2000, 100, 1000, 10000,
                                                      true);
  fcl::broadphaseCollisionTestBinaryAABB_Tree<double>(2000, 1000, 1000, 10000,
                                                      true);
#else
  fcl::broadphaseCollisionTestBinaryAABB_Tree<double>(2000, 10, 100, 1, false);
  fcl::broadphaseCollisionTestBinaryAABB_Tree<double>(2000, 100, 100, 1, false);
#endif
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
