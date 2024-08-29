#include <gtest/gtest.h>

#include "fcl/narrowphase/colliding_volume.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

namespace fcl {

constexpr int test_size = 10000;
constexpr double max_volume = std::numeric_limits<double>::max();
constexpr std::size_t max_n_contact = std::numeric_limits<std::size_t>::max();

// Calculate the volume of the intersection of two spheres
template <typename S>
S calcIntersectionVolumeForTwoSpheres(S radius1, S radius2,
                                      const Vector3<S>& translation) {
  if (radius1 > radius2) std::swap(radius1, radius2);
  const S distance = translation.norm();
  if (distance >= (radius2 + radius1)) return 0;
  if (distance <= (radius2 - radius1))
    return (S)4.0 * constants<S>::pi() * radius1 * radius1 * radius1 / (S)3.0;

  const S cos_1 =
      (radius1 * radius1 + distance * distance - radius2 * radius2) /
      (2 * radius1 * distance);
  const S cos_2 =
      (radius2 * radius2 + distance * distance - radius1 * radius1) /
      (2 * radius2 * distance);
  const S h_1 = radius1 * (1 - cos_1);
  const S h_2 = radius2 * (1 - cos_2);

  const auto pi_3 = constants<S>::pi() / 3;
  const S volume = pi_3 * ((3 * radius2 - h_2) * h_2 * h_2 +
                           (3 * radius1 - h_1) * h_1 * h_1);
  return volume;
}

// Calculate the area where two circles intersect
template <typename S>
S calcIntersectionAreaForTwoCircles(const S radius1, const S radius2,
                                    const S distance) {
  if (distance > radius2 + radius1) return 0;
  if (distance < abs(radius2 - radius1))
    return radius1 > radius2 ? radius2 * radius2 * constants<S>::pi()
                             : radius1 * radius1 * constants<S>::pi();

  const S theta_1 =
      acos((radius1 * radius1 + distance * distance - radius2 * radius2) /
           (2 * radius1 * distance));
  const S theta_2 =
      acos((radius2 * radius2 + distance * distance - radius1 * radius1) /
           (2 * radius2 * distance));
  const S area1 = radius1 * radius1 * (theta_1 - sin(theta_1) * cos(theta_1));
  const S area2 = radius2 * radius2 * (theta_2 - sin(theta_2) * cos(theta_2));
  return area1 + area2;
}

template <typename S>
S cylinderCollisionVolumeLowerBound(
    const std::shared_ptr<Cylinder<S>>& cylinder1,
    const std::shared_ptr<Cylinder<S>>& cylinder2,
    const Vector3<S>& translation) {
  const S height =
      std::min(cylinder1->lz / 2, cylinder2->lz / 2 + translation.z()) -
      std::max(-cylinder1->lz / 2, -cylinder2->lz / 2 + translation.z());
  if (height < 0) return 0;
  const S area = calcIntersectionAreaForTwoCircles(
      cylinder1->radius, cylinder2->radius,
      sqrt(translation.x() * translation.x() +
           translation.y() * translation.y()));
  return height * area;
}

template <typename S>
void triangleAABBCollisionVolumeTest() {
  for (int j = 0; j < test_size; j++) {
    constexpr int vertices_size = 3;
    std::vector<Vector3<S>> vertices;
    for (int i = 0; i < vertices_size; i++) {
      vertices.push_back(Vector3<S>::Random());
    }

    const S bvhTriangleThickness = 0.001;
    AABB<S> aabb_output;
    Transform3<S> tf_AABB_output;
    detail::calcTriangleAABB<S>(vertices[0], vertices[1], vertices[2],
                                bvhTriangleThickness, aabb_output,
                                tf_AABB_output);
    const S tolerance = 1e-6;
    const Vector3<S> delta{tolerance, tolerance, tolerance};
    aabb_output.expand(delta);

    const auto point_to_origin = tf_AABB_output.inverse();

    Vector3<S> normal =
        (vertices[0] - vertices[1]).cross(vertices[0] - vertices[2]);
    normal.normalize();
    for (int i = 0; i < vertices_size; i++) {
      EXPECT_TRUE(aabb_output.contain(point_to_origin * vertices[0]));
      EXPECT_TRUE(aabb_output.contain(
          point_to_origin *
          (vertices[i] + normal * (0.5 * bvhTriangleThickness - tolerance))));
      EXPECT_TRUE(aabb_output.contain(
          point_to_origin *
          (vertices[i] - normal * (0.5 * bvhTriangleThickness - tolerance))));
      EXPECT_FALSE(aabb_output.contain(
          point_to_origin *
          (vertices[i] + (0.5 + 1e-2) * normal * bvhTriangleThickness)));
      EXPECT_FALSE(aabb_output.contain(
          point_to_origin *
          (vertices[i] - (0.5 + 1e-2) * normal * bvhTriangleThickness)));
    }

    auto testTriangleSide = [&aabb_output, &point_to_origin](
                                const Vector3<S>& vertice_lhi,
                                const Vector3<S>& vertice_rhi) {
      EXPECT_TRUE(aabb_output.contain(point_to_origin *
                                      (0.5 * vertice_lhi + 0.5 * vertice_rhi)));
      EXPECT_FALSE(aabb_output.contain(
          point_to_origin * ((1 + 1e-3) * vertice_rhi - 1e-3 * vertice_lhi)));
      EXPECT_FALSE(aabb_output.contain(
          point_to_origin * ((1 + 1e-3) * vertice_lhi - 1e-3 * vertice_rhi)));
    };

    testTriangleSide(vertices[0], vertices[1]);
    testTriangleSide(vertices[1], vertices[2]);
    testTriangleSide(vertices[2], vertices[0]);
  }
}

template <typename S>
void axisAlignedBoxPairCollisionVolumeTest() {
  std::vector<Vector3<S>> sides{
      Vector3<S>(1, 2, 3),
      Vector3<S>(1, 1, 1),
      Vector3<S>(2, 1, 2),
      Vector3<S>(2, 1, 3),
  };

  const std::array<S, 6> extents{-1, -1, -1, 1, 1, 1};
  for (const auto& side1 : sides) {
    for (const auto& side2 : sides) {
      auto box1 = std::make_shared<Box<S>>(side1);
      box1->computeLocalAABB();
      auto box2 = std::make_shared<Box<S>>(side2);
      box2->computeLocalAABB();
      for (int i = 0; i < test_size; i++) {
        Transform3<S> tf1;
        test::generateRandomTransform(extents, tf1);
        Transform3<S> tf2;
        tf2.setIdentity();
        Vector3<S> box1_to_box2_translation = Vector3<S>::Random();
        tf2.translation() = box1_to_box2_translation;
        tf2 = tf1 * tf2;

        CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
        CollidingVolumeResult<S> volume_result;
        approximateCollisionVolume(box1.get(), tf1, box2.get(), tf2,
                                   volume_request, volume_result);

        const AABB<S> aabb_after_transform =
            translate(box2->aabb_local, box1_to_box2_translation);
        AABB<S> overlap_part;
        double volume = 0;
        if (box1->aabb_local.overlap(aabb_after_transform, overlap_part))
          volume = overlap_part.volume();
        EXPECT_TRUE(abs(volume_result.totalVolume() - volume) < 1e-5);
      }
    }
  }
}

template <typename S>
void sphereLowerBoundCollisionVolumeTest() {
  std::vector<S> radius{1, 1.1, 1.3, 1.8, 2};
  const std::array<S, 6> extents{-1, -1, -1, 1, 1, 1};
  for (const auto r1 : radius) {
    for (const auto r2 : radius) {
      auto sphere1 = std::make_shared<Sphere<S>>(r1);
      sphere1->computeLocalAABB();
      auto sphere2 = std::make_shared<Sphere<S>>(r2);
      sphere2->computeLocalAABB();
      for (int i = 0; i < test_size; i++) {
        Transform3<S> tf1;
        test::generateRandomTransform(extents, tf1);
        Transform3<S> tf2;
        test::generateRandomTransform(extents, tf2);

        CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
        CollidingVolumeResult<S> volume_result;
        approximateCollisionVolume(sphere1.get(), tf1, sphere2.get(), tf2,
                                   volume_request, volume_result);
        const Vector3<S> translation = tf1.translation() - tf2.translation();
        const auto lowerBound = calcIntersectionVolumeForTwoSpheres(
            sphere1->radius, sphere2->radius, translation);
        EXPECT_TRUE(volume_result.totalVolume() >= lowerBound);
      }
    }
  }
}

template <typename S>
void cylinderLowerBoundCollisionVolumeTest() {
  auto cylinder1 = std::make_shared<Cylinder<S>>(1.6, 1.0);
  cylinder1->computeLocalAABB();
  auto cylinder2 = std::make_shared<Cylinder<S>>(1.5, 1.0);
  cylinder2->computeLocalAABB();

  const std::array<S, 6> extents{-1, -1, -1, 1, 1, 1};
  for (int i = 0; i < test_size; i++) {
    Transform3<S> tf1;
    test::generateRandomTransform(extents, tf1);
    Transform3<S> tf2;
    tf2.setIdentity();
    Vector3<S> cylinder2_to_cylinder2_translation = Vector3<S>::Random();
    tf2.translation() = cylinder2_to_cylinder2_translation;
    tf2 = tf1 * tf2;

    CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
    CollidingVolumeResult<S> volume_result;
    approximateCollisionVolume(cylinder1.get(), tf1, cylinder2.get(), tf2,
                               volume_request, volume_result);
    const auto lowerBound = cylinderCollisionVolumeLowerBound(
        cylinder1, cylinder2, cylinder2_to_cylinder2_translation);
    if (volume_result.totalVolume() < lowerBound) {
      std::cout << "tf1" << std::endl;
      std::cout << tf1.matrix() << std::endl;
      std::cout << "tf2" << std::endl;
      std::cout << tf2.matrix() << std::endl;
      std::cout << "*****************" << std::endl;
    }
    EXPECT_TRUE(volume_result.totalVolume() >= lowerBound);
  }
}

template <typename Shape1, typename Shape2>
void shapePairVolumeInstanceTest(const Shape1& shape1, const Shape2& shape2) {
  static_assert(std::is_same<typename Shape1::S, typename Shape2::S>::value,
                "Test shape mismatch");
  using S = typename Shape1::S;

  fcl::Transform3<S> shape1_pose;
  {
    shape1_pose.setIdentity();
    shape1_pose.translation().x() = 0.530871756374835968;
    shape1_pose.translation().y() = 0.361621372401714325;
    shape1_pose.translation().z() = -0.2295833881944417953;

    // The rotation matrix
    fcl::Matrix3<S> rotation_matrix;
    rotation_matrix.row(0) = fcl::Vector3<S>(
        0.1465600677548429542, -0.4152570123227734, 0.8978205612796243962);
    rotation_matrix.row(1) = fcl::Vector3<S>(
        0.9092629872412127945, 0.4139898295249834215, 0.04304928667308412921);
    rotation_matrix.row(2) = fcl::Vector3<S>(
        -0.3895650992746680918, 0.8100456991840886412, 0.4382522089624635298);
    shape1_pose.linear().matrix() = rotation_matrix;
  }

  fcl::Transform3<S> shape2_pose = shape1_pose;
  {
    shape2_pose.translation().x() = 0.2356679795984277237;
    shape2_pose.translation().y() = -0.05616388602768807026;
    shape2_pose.translation().z() = 0.8171584263755837796;
  }

  CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
  CollidingVolumeResult<S> volume_result;
  approximateCollisionVolume(&shape1, shape1_pose, &shape2, shape2_pose,
                             volume_request, volume_result);
  EXPECT_TRUE(volume_result.numContacts() > 0);
}

template <typename S>
void cylinderBoxLowerBoundCollisionVolumeTest() {
  srand(0);
  const S length = 1.0;
  const S width = 0.1 * length;
  const S height = 1.0;
  auto cylinder = std::make_shared<Cylinder<S>>(length, height);
  cylinder->computeLocalAABB();

  auto box = std::make_shared<Box<S>>(length, width, height);
  box->computeLocalAABB();
  for (int i = 0; i < test_size; i++) {
    Transform3<S> tf;
    tf.setIdentity();
    tf.translation() =
        Vector3<S>(test::rand_interval<S>(0.5 * length, 1.5 * length), 0, 0);

    CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
    CollidingVolumeResult<S> volume_result;
    approximateCollisionVolume(cylinder.get(), Transform3<S>::Identity(),
                               box.get(), tf, volume_request, volume_result);
    const auto lowerBound =
        width * height * (1.5 * length - tf.translation().x());
    EXPECT_TRUE((volume_result.totalVolume() + 1e-6) > lowerBound);
  }
}

template <typename S>
void coneBoxLowerBoundCollisionVolumeTest() {
  const S length = 1.0;
  const S width = 0.1 * length;
  const S height = 1.0;
  std::shared_ptr<Cone<S>> cone = std::make_shared<Cone<S>>(length, height);
  cone->computeLocalAABB();

  std::shared_ptr<Box<S>> box = std::make_shared<Box<S>>(length, width, height);
  box->computeLocalAABB();
  for (int i = 0; i < test_size; i++) {
    Transform3<S> tf;
    tf.setIdentity();
    tf.translation() =
        Vector3<S>(test::rand_interval<S>(0.5 * length, 1.5 * length), 0, 0);

    CollidingVolumeRequest<S> volume_request(max_n_contact, max_volume);
    CollidingVolumeResult<S> volume_result;
    approximateCollisionVolume(cone.get(), Transform3<S>::Identity(), box.get(),
                               tf, volume_request, volume_result);

    const S bottom_length = (1.5 * length - tf.translation().x());
    const S height_length = bottom_length / length * height;
    const auto lowerBound = bottom_length * height_length * width / 2;
    EXPECT_TRUE(volume_result.totalVolume() > lowerBound);
  }
}

}  // namespace fcl

GTEST_TEST(CollisionVolumeTest, AxisAlignedBoxPairTest) {
  fcl::axisAlignedBoxPairCollisionVolumeTest<float>();
  fcl::axisAlignedBoxPairCollisionVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, SphereLowerBoundTest) {
  fcl::sphereLowerBoundCollisionVolumeTest<float>();
  fcl::sphereLowerBoundCollisionVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, CylinderLowerBoundTest) {
  /*
  using S = double;
  auto cylinder1 = std::make_shared<fcl::Cylinder<S>>(1.6, 1.0);
  cylinder1->computeLocalAABB();
  auto cylinder2 = std::make_shared<fcl::Cylinder<S>>(1.5, 1.0);
  cylinder2->computeLocalAABB();
  fcl::shapePairVolumeInstanceTest(*cylinder1, *cylinder2);*/

  fcl::cylinderLowerBoundCollisionVolumeTest<float>();
  fcl::cylinderLowerBoundCollisionVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, CylinderBoxLowerBoundTest) {
  fcl::cylinderBoxLowerBoundCollisionVolumeTest<float>();
  fcl::cylinderBoxLowerBoundCollisionVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, ConeBoxLowerBoundTest) {
  fcl::coneBoxLowerBoundCollisionVolumeTest<float>();
  fcl::coneBoxLowerBoundCollisionVolumeTest<double>();
}

GTEST_TEST(CollisionVolumeTest, TriangleAABBTest) {
  fcl::triangleAABBCollisionVolumeTest<float>();
  fcl::triangleAABBCollisionVolumeTest<double>();
}

int main(int argc, char* argv[]) {
  std::cout.precision(std::numeric_limits<float>::max_digits10 + 10);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
