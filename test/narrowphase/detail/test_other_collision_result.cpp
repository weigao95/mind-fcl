//
// Created by wei on 23-3-16.
//

#include <gtest/gtest.h>

#include "fcl/geometry/geometric_shape_to_BVH_model.h"
#include "fcl/narrowphase/collision.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

namespace fcl {

template <typename S>
void simpleContactResultTest() {
  srand(2);
  std::vector<CollisionObject<S>*> env;
  test::generateEnvironments<S>(env, 200.0, 4);

  // TODO(wei): revive this test with octree2
  /*const auto octree =
      std::shared_ptr<const octomap::OcTree>(test::generateOcTree(1.0));
  auto tree_geometry = std::make_shared<OcTree<S>>(octree);
  tree_geometry->computeLocalAABB();

  Transform3<S> tf_1;
  tf_1.setIdentity();
  Transform3<S> tf_2 = tf_1;
  tf_2.translation() = Vector3<S>(0, 0, 0.5);
  for (std::size_t request_count : {1, 10, 100, 1000}) {
    fcl::CollisionRequest<S> request(request_count);

    // Two different result
    fcl::CollisionResult<S> default_result;
    fcl::collide(tree_geometry.get(), tf_1, tree_geometry.get(), tf_2, request,
                 default_result);

    std::size_t counter = 0;
    auto counting_processor = [&counter, request_count](
                                  const Contact<S>&, bool& keep_this,
                                  bool& user_stop) -> void {
      keep_this = false;
      counter += 1;
      user_stop = (counter >= request_count);
    };

    CollisionResult<S> counting_result(std::move(counting_processor));
    fcl::collide(tree_geometry.get(), tf_1, tree_geometry.get(), tf_2, request,
                 counting_result);
    EXPECT_EQ(default_result.numContacts(), counter);
    EXPECT_EQ(counting_result.numContacts(), 0u);
  }*/
}

}  // namespace fcl

GTEST_TEST(CollisionResultTest, CountingContactTest) {
  fcl::simpleContactResultTest<float>();
  fcl::simpleContactResultTest<double>();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}