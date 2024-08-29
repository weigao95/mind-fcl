//
// Created by Wei Gao on 2024/8/8.
//

#pragma once

#include <cassert>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <vector>

namespace fcl {
namespace detail {

template <typename Object>
class SimpleVectorObjectAllocator {
 public:
  static_assert(std::is_default_constructible<Object>::value,
                "Object must be default constructible");
  static_assert(std::is_trivially_destructible<Object>::value,
                "Object must be trivially destructible");
  using Index = std::uint32_t;
  static constexpr Index kInvalidIndex = 0xffffffff;

  // Allocate an object
  Index AllocateObject();
  template <typename... Parameters>
  Index ConstructObject(Parameters&&... parameters);

  // Release an object for further allocation
  void Free(Index index);

  // Query constructed objects
  Object& Get(Index index) { return buffer_[index]; }
  const Object& Get(Index index) const { return buffer_[index]; }

  // Query the vector
  const std::vector<Object>& raw_vector() const { return buffer_; }
 private:
  std::vector<Object> buffer_;

  // For free
  Index free_begin_{kInvalidIndex};
  std::vector<Index> free_index_list_;
};

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2 = 10>
class FixedSizeFreeList {
 private:
  struct ObjectStorage {
    // Object data
    Object object;

    // If this is an allocated (and valid) object, then this is the index of
    // the object.
    // else (if the object is freed), this is the index of the next freed
    // object (which forms a linked list)
    std::uint32_t index;
  };

 public:
  using Index = std::uint32_t;
  static constexpr Index kInvalidIndex = 0xffffffff;
  static constexpr Index kPageSize = 1 << kPageSizeAsPowerOf_2;
  static constexpr Index kPageSizeInObjects = kPageSize;
  static constexpr Index kPageOffsetShift = kPageSizeAsPowerOf_2;
  static constexpr Index kInPageIndexMask = kPageSize - 1;
  static_assert(std::is_default_constructible<Object>::value,
                "Object must be default constructible");
  static_assert(std::is_trivially_destructible<Object>::value,
                "Object must be trivially destructible");

  // Construction
  explicit FixedSizeFreeList(Index max_num_pages);
  FixedSizeFreeList(const FixedSizeFreeList&) = delete;
  ~FixedSizeFreeList();

  // Insert object and return a handle of that object
  template <typename... Parameters>
  Index ConstructObject(Parameters&&... parameters);
  Index AllocateObject();

  // Release an object for further allocation
  void Free(Index index);

  // Query constructed objects
  Object& Get(Index index) { return getStorage(index).object; }
  const Object& Get(Index index) const { return getStorage(index).object; }

  // clang-format off
  inline Index capacity() const { return kPageSizeInObjects * page_table_.max_num_pages; }
  inline Index capacity_in_objects() const { return capacity(); }
  // clang-format on

 private:
  // Size info
  struct {
    Index max_num_pages{0};
    Index num_allocated_objects{0};
    Index num_allocated_pages{0};
    std::vector<std::unique_ptr<ObjectStorage[]>> pages{};
    Index freed_object_list_{kInvalidIndex};
  } page_table_;

  // clang-format off
  void initializePageTable(Index max_num_pages);
  ObjectStorage& getStorage(Index index) { return page_table_.pages[index >> kPageOffsetShift].get()[index & kInPageIndexMask]; };
  const ObjectStorage& getStorage(Index index) const { return page_table_.pages[index >> kPageOffsetShift].get()[index & kInPageIndexMask]; }
  // clang-format on
};

}  // namespace detail
}  // namespace fcl

#include "fcl/broadphase/binary_AABB_tree_allocator-inl.h"
