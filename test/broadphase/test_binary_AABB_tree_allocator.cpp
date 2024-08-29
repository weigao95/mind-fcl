//
// Created by Wei Gao on 2024/8/8.
//
#include <gtest/gtest.h>

#include "fcl/broadphase/binary_AABB_tree_allocator.h"

namespace fcl {
namespace detail {

template <typename UintAllocator>
void testUintAllocator(UintAllocator& allocator, std::uint32_t test_n,
                       bool invalid_beyond_n) {
  for (std::uint32_t i = 0; i < test_n; i++) {
    allocator.ConstructObject(i * i + 1);
  }

  if (invalid_beyond_n) {
    for (std::uint32_t i = 0; i < test_n; i++) {
      auto output = allocator.ConstructObject(i * i + 2);
      EXPECT_EQ(output, 0xffffffff);
    }
  }

  for (std::uint32_t i = 0; i < test_n; i++) {
    auto elem = allocator.Get(i);
    EXPECT_EQ(elem, i * i + 1);
  }

  for (std::uint32_t i = 0; i < test_n; i++) {
    allocator.Free(i);
  }

  for (std::uint32_t i = 0; i < test_n; i++) {
    auto output = allocator.ConstructObject(i * i + 2);
    EXPECT_NE(output, 0xffffffff);
    auto elem = allocator.Get(output);
    EXPECT_EQ(elem, i * i + 2);
  }
}

void testSimpleVectorAllocator(std::uint32_t test_n) {
  SimpleVectorObjectAllocator<std::uint32_t> allocator;
  testUintAllocator(allocator, test_n, false);
  const auto& raw_vec = allocator.raw_vector();
  EXPECT_EQ(raw_vec.size(), test_n);
}

template <std::uint32_t kPageSizeAsPowerOf_2>
void testFixedSizeFreeListSimple(std::uint32_t max_n_pages) {
  FixedSizeFreeList<std::uint32_t, kPageSizeAsPowerOf_2> allocator(max_n_pages);
  EXPECT_EQ(allocator.capacity_in_objects(),
            max_n_pages * (1 << kPageSizeAsPowerOf_2));
  testUintAllocator(allocator, allocator.capacity(), true);
}

}  // namespace detail
}  // namespace fcl

GTEST_TEST(FixedSizeFreeListTest, SimpleVectorTest) {
  fcl::detail::testSimpleVectorAllocator(1);
  fcl::detail::testSimpleVectorAllocator(2);
  fcl::detail::testSimpleVectorAllocator(1024);
  fcl::detail::testSimpleVectorAllocator(1024 * 1024);
}

GTEST_TEST(FixedSizeFreeListTest, SimpleTest) {
  fcl::detail::testFixedSizeFreeListSimple<1>(1);
  fcl::detail::testFixedSizeFreeListSimple<1>(2);
  fcl::detail::testFixedSizeFreeListSimple<10>(1024);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}