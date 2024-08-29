#pragma once

namespace fcl {
namespace detail {

template <typename Object>
std::uint32_t SimpleVectorObjectAllocator<Object>::AllocateObject() {
  if (free_begin_ != kInvalidIndex) {
    const auto index = free_begin_;
    assert(index < free_index_list_.size());
    free_begin_ = free_index_list_[index];
    return index;
  }

  // Need to allocate a new one
  const auto index = static_cast<Index>(buffer_.size());
  buffer_.emplace_back();
  free_index_list_.emplace_back(index);
  return index;
}

template <typename Object>
template <typename... Parameters>
std::uint32_t SimpleVectorObjectAllocator<Object>::ConstructObject(
    Parameters&&... parameters) {
  const auto index = AllocateObject();
  if (index == kInvalidIndex) return index;

  // Construct the element
  auto& storage = Get(index);
  ::new (&storage) Object(std::forward<Parameters>(parameters)...);
  return index;
}

template <typename Object>
void SimpleVectorObjectAllocator<Object>::Free(Index index) {
  assert(index < buffer_.size());
  assert(buffer_.size() == free_index_list_.size());
  free_index_list_[index] = free_begin_;
  free_begin_ = index;
}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::FixedSizeFreeList(
    Index max_num_pages) {
  static_assert(kPageSizeAsPowerOf_2 >= 1, "Invalid page size");
  initializePageTable(max_num_pages);
}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::~FixedSizeFreeList() {}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
void FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::initializePageTable(
    Index max_num_pages) {
  // Does not allow double init
  if (page_table_.max_num_pages > 0) {
    return;
  }

  // Set the meta and allocate space
  page_table_.max_num_pages = max_num_pages;
  page_table_.num_allocated_objects = 0;
  page_table_.num_allocated_pages = 0;

  // Allocate the page table
  page_table_.pages.resize(max_num_pages);
  for (Index i = 0; i < max_num_pages; i++) page_table_.pages[i] = nullptr;
}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
void FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::Free(Index index) {
  auto& storage = getStorage(index);
  // Not invoke destructor as it is default-destruct

  storage.index = page_table_.freed_object_list_;
  page_table_.freed_object_list_ = index;
}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
std::uint32_t
FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::AllocateObject() {
  if (page_table_.freed_object_list_ != kInvalidIndex) {
    // Get the location
    const auto construct_at = page_table_.freed_object_list_;
    auto& storage = getStorage(construct_at);

    // Update the freed list
    const auto next_free = storage.index;
    page_table_.freed_object_list_ = next_free;

    // Use this as the new elem
    storage.index = construct_at;
    return construct_at;
  }

  // Request a new location
  const auto construct_at_offset = page_table_.num_allocated_objects;
  if (construct_at_offset >= page_table_.num_allocated_pages * kPageSize) {
    // Check capacity
    const auto new_page_index = page_table_.num_allocated_pages;
    if (new_page_index >= page_table_.max_num_pages) {
      return kInvalidIndex;
    }

    // Allocate a new page
    page_table_.pages[new_page_index] =
        std::unique_ptr<ObjectStorage[]>(new ObjectStorage[kPageSize]);

    // Update the page count
    page_table_.num_allocated_pages += 1;
  }

  // Make index
  const auto construct_at = construct_at_offset;
  assert((construct_at >> kPageOffsetShift) + 1 ==
         page_table_.num_allocated_pages);

  // Construct the new elem
  auto& storage = getStorage(construct_at);
  storage.index = construct_at;

  // Update the index and done
  page_table_.num_allocated_objects += 1;
  return construct_at;
}

template <typename Object, std::uint32_t kPageSizeAsPowerOf_2>
template <typename... Parameters>
std::uint32_t FixedSizeFreeList<Object, kPageSizeAsPowerOf_2>::ConstructObject(
    Parameters&&... parameters) {
  const auto index = AllocateObject();
  if (index == kInvalidIndex) return index;

  // Construct the element
  auto& storage = getStorage(index);
  ::new (&storage.object) Object(std::forward<Parameters>(parameters)...);
  assert(storage.index == index);
  return index;
}

}  // namespace detail
}  // namespace fcl
