//
// Created by Wei Gao on 2024/8/9.
//

#pragma once

#include <stack>

#include "fcl/broadphase/broadphase_common.h"
#include "fcl/broadphase/binary_AABB_tree_allocator.h"

namespace fcl {
namespace detail {

template <typename S, template <typename Object> class ObjectAllocator>
class BinaryAABB_Tree {
  // Copy the allocator info
  using AllocatorIndex = std::uint32_t;
  static constexpr AllocatorIndex kInvalidAllocatorIndex = 0xffffffff;

  // Typedef for node status
  struct NodeStatus {
    // Inner or leaf
    // clang-format off
    static constexpr std::uint32_t kLeafOrInnerBit = (1 << 0);
    static constexpr std::uint32_t kLeafBitValue = 0;
    inline bool IsLeaf() const { return (status & kLeafOrInnerBit) == kLeafBitValue; };
    inline bool IsInner() const { return (status & kLeafOrInnerBit) != kLeafBitValue; }
    inline void SetAsLeaf() { status &= (~kLeafOrInnerBit); };
    inline void SetAsInner() { status |= kLeafOrInnerBit; }

    // Removed
    static constexpr std::uint32_t kRemovedBit = (1 << 1);
    static constexpr std::uint32_t kNotRemovedBitValue = 0;
    inline bool IsRemoved() const { return (status & kRemovedBit) != kNotRemovedBitValue; }
    inline void SetAsRemoved() { status |= kRemovedBit; }
    // clang-format on

    // Actual status
    static constexpr std::uint32_t kDefaultStatus = 0;
    std::uint32_t status{kDefaultStatus};
  };

  // Typedef for node
  struct Node {
    // BV of the node
    AABB<S> bv{};

    // Children of user id
    union {
      std::array<AllocatorIndex, 2> children;
      std::uint64_t user_id;
    };

    // Parent and meta
    AllocatorIndex parent{kInvalidAllocatorIndex};
    NodeStatus status{};
  };

  // Typedef for allocator
  using Allocator = ObjectAllocator<Node>;
  Allocator node_allocator_;

  // Root node
  AllocatorIndex root_node_{kInvalidAllocatorIndex};

  // Mapping from user_id to leaf node index
  using UserIdMap = std::unordered_map<std::uint64_t, AllocatorIndex>;
  UserIdMap user_id_map_;

 public:
  explicit BinaryAABB_Tree();
  ~BinaryAABB_Tree() = default;

  // Build the tree, which might use the internal node allocator of an
  // existing AABB tree. In that case, the BuildTreeExternal thread can run
  // concurrently with the collision detection methods (Visit,
  // SingleObject/Tree/Self-Collision) of the existing AABB tree, despite
  // the read/write access. However, multiple BuildTreeExternal threads
  // cannot run concurrently.
  // The objects would be MUTATED, as it is both input and algorithm BUFFER
  static AllocatorIndex BuildTreeExternal(BroadphaseObjectInfo<S>* objects,
                                          std::uint32_t n_objects,
                                          Allocator& allocator,
                                          UserIdMap* user_id_map = nullptr);
  static std::uint32_t DestructTreeExternal(AllocatorIndex root_node,
                                            Allocator& allocator);

  // Build the tree in place, which internally invoked the BuildTreeExternal
  // and DestructTreeExternal. This method can NOT be used concurrently with
  // collision detection methods below.
  void Rebuild(BroadphaseObjectInfo<S>* objects, std::uint32_t n_objects);

  // For updating the tree structure, which invoke BuildTreeExternal and
  // DestructTreeExternal internally.
  // As mentioned above, PrepareUpdateStructure can be invoked concurrently
  // with visit and collision detection methods in aother ONE thread. However, 
  // Apply/Abort-Update can NOT be invoked concurrently with collision 
  // detection methods.
  struct TreeUpdateState {
    AllocatorIndex root_index{kInvalidAllocatorIndex};
    UserIdMap user_id_map;
  };
  void PrepareUpdateStructure(TreeUpdateState& state);
  void ApplyUpdateStructure(TreeUpdateState&& state);
  void AbortUpdateStructure(TreeUpdateState& state);
  void PrepareAddNewObjects(BroadphaseObjectInfo<S>* new_objects,
                            std::uint32_t n_new_objects,
                            TreeUpdateState& state);
  bool ApplyAddNewObjects(TreeUpdateState& state);
  void AbortAddNewObjects(TreeUpdateState& state);

  // Object-wise update method
  bool RemoveObject(std::uint64_t object_user_id);
  bool UpdateObjectAABB(std::uint64_t object_user_id, const AABB<S>& new_AABB);

  // Visit the tree with a given functor, return 1) whether current node (and
  // all its children can be terminated or not); 2) Overall termination
  using VisitorFn = std::function<void(
      const AABB<S>& node_aabb, bool is_leaf, std::uint64_t user_id_if_leaf,
      bool& current_node_done, bool& overall_done)>;
  void VisitTree(const VisitorFn& visitor_fn) const;

  // Collision function with another object
  using CollisionFn =
      std::function<bool(std::uint64_t leaf1_user_id,
                         std::uint64_t leaf2_user_id, void* collision_fn_data)>;
  void SingleObjectCollision(const BroadphaseObjectInfo<S>& object,
                             const CollisionFn& collision_fn,
                             void* collision_fn_data) const;
  void TreeCollision(const BinaryAABB_Tree<S, ObjectAllocator>& tree2,
                     const CollisionFn& collision_fn,
                     void* collision_fn_data) const;
  void SelfCollision(const CollisionFn& collision_fn,
                     void* collision_fn_data) const;

  // State query
  Allocator& allocator() { return node_allocator_; };
  std::size_t n_leaves() const { return user_id_map_.size(); }

  // State checking
  bool SanityCheck() const;
};

}  // namespace detail
}  // namespace fcl

#include "fcl/broadphase/binary_AABB_tree-inl.h"
