#pragma once

namespace fcl {
namespace detail {

template <typename S>
bool sortBroadphaseObjectByCenterCompare(const BroadphaseObjectInfo<S>& a,
                                         const BroadphaseObjectInfo<S>& b,
                                         int axis) {
  return a.bv.center()[axis] < b.bv.center()[axis];
}

template <typename S, template <typename Object> class ObjectAllocator>
BinaryAABB_Tree<S, ObjectAllocator>::BinaryAABB_Tree()
    : node_allocator_(), root_node_(kInvalidAllocatorIndex) {
  assert(user_id_map_.empty());
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::Rebuild(
    BroadphaseObjectInfo<S>* objects, std::uint32_t n_objects) {
  DestructTreeExternal(root_node_, node_allocator_);
  root_node_ =
      BuildTreeExternal(objects, n_objects, node_allocator_, &user_id_map_);
}

template <typename S, template <typename Object> class ObjectAllocator>
std::uint32_t BinaryAABB_Tree<S, ObjectAllocator>::DestructTreeExternal(
    AllocatorIndex root_node, Allocator& allocator) {
  // Special case: nothing to remove
  if (root_node == kInvalidAllocatorIndex) {
    return 0;
  }

  // Processing loop
  std::uint32_t n_removed = 0U;
  std::stack<AllocatorIndex> task_stack;
  task_stack.push(root_node);
  while (!task_stack.empty()) {
    const auto node_index = task_stack.top();
    task_stack.pop();

    // Obtain the left and right
    assert(node_index != kInvalidAllocatorIndex);
    AllocatorIndex left = kInvalidAllocatorIndex;
    AllocatorIndex right = kInvalidAllocatorIndex;
    {
      const auto& node = allocator.Get(node_index);
      if (node.status.IsInner()) {
        left = node.children[0];
        right = node.children[1];
      }
    }

    // Release the current
    allocator.Free(node_index);
    n_removed += 1U;
    if (left != kInvalidAllocatorIndex) {
      assert(right != kInvalidAllocatorIndex);
      task_stack.push(left);
      task_stack.push(right);
    }
  }

  // Done
  return n_removed;
}

template <typename S, template <typename Object> class ObjectAllocator>
std::uint32_t BinaryAABB_Tree<S, ObjectAllocator>::BuildTreeExternal(
    BroadphaseObjectInfo<S>* objects, std::uint32_t n_objects,
    Allocator& allocator, UserIdMap* user_id_map) {
  // Special case: no objects input
  if (n_objects == 0) {
    return kInvalidAllocatorIndex;
  }

  // Clear existing map if present
  if (user_id_map != nullptr) {
    user_id_map->clear();
  }

  // General case
  struct StackElement {
    std::uint32_t begin;
    std::uint32_t n_object;
    AllocatorIndex parent;
    bool is_left_child_of_parent;
  };
  std::stack<StackElement> task_stack;
  task_stack.push({0, n_objects, kInvalidAllocatorIndex, false});

  // Processing loop
  AllocatorIndex root_index{kInvalidAllocatorIndex};
  while (!task_stack.empty()) {
    // Current element
    const auto task = task_stack.top();
    task_stack.pop();

    // Leaf case of only one objects
    if (task.n_object == 1U) {
      // Make a new node
      const auto node_index = allocator.AllocateObject();
      if (node_index == kInvalidAllocatorIndex) {
        return node_index;
      }

      // Assign meta
      auto& node = allocator.Get(node_index);
      const auto& object = objects[task.begin];
      node.bv = object.bv;
      node.user_id = object.user_id;
      node.parent = task.parent;
      node.status.SetAsLeaf();

      // Insert to user id map
      if (user_id_map != nullptr) {
        user_id_map->insert(std::make_pair(node.user_id, node_index));
      }

      // Assign meta of parent
      if (task.parent != kInvalidAllocatorIndex) {
        auto& parent_node = allocator.Get(task.parent);
        if (task.is_left_child_of_parent) {
          parent_node.children[0] = node_index;
        } else {
          parent_node.children[1] = node_index;
        }
      } else {
        assert(root_index == kInvalidAllocatorIndex);  // Only once
        root_index = node_index;
      }

      // Done with this node
      continue;
    }

    // Else
    assert(task.n_object >= 2);
    AABB<S> node_bv = objects[task.begin].bv;
    for (std::uint32_t i = 1; i < task.n_object; i++) {
      const AABB<S>& bv_i = objects[task.begin + i].bv;
      node_bv += bv_i;
    }

    // Select the best axis
    const auto n_center = task.n_object / 2;
    {
      int split_axis = 0;
      S extent[3] = {node_bv.width(), node_bv.height(), node_bv.depth()};
      if (extent[1] > extent[0]) split_axis = 1;
      if (extent[2] > extent[split_axis]) split_axis = 2;

      assert(n_center >= 1);
      std::nth_element(objects + task.begin, objects + (task.begin + n_center),
                       objects + (task.begin + task.n_object),
                       std::bind(&sortBroadphaseObjectByCenterCompare<S>,
                                 std::placeholders::_1, std::placeholders::_2,
                                 std::ref(split_axis)));
    }

    // Make node
    const auto node_index = allocator.AllocateObject();
    if (node_index == kInvalidAllocatorIndex) {
      return node_index;
    }

    // Assign meta
    auto& new_node = allocator.Get(node_index);
    new_node.bv = node_bv;
    new_node.parent = task.parent;
    new_node.status.SetAsInner();

    // Assign meta of parent
    if (task.parent != kInvalidAllocatorIndex) {
      auto& parent_node = allocator.Get(task.parent);
      if (task.is_left_child_of_parent) {
        parent_node.children[0] = node_index;
      } else {
        parent_node.children[1] = node_index;
      }
    } else {
      assert(root_index == kInvalidAllocatorIndex);  // Only once
      root_index = node_index;
    }

    // Push new task
    task_stack.push(
        {task.begin + n_center, task.n_object - n_center, node_index, false});
    task_stack.push({task.begin, n_center, node_index, true});
  }

  // Done
  return root_index;
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::PrepareUpdateStructure(
    TreeUpdateState& state) {
  // Special case of empty
  if (root_node_ == kInvalidAllocatorIndex) {
    state.root_index = kInvalidAllocatorIndex;
    state.user_id_map.clear();
    return;
  }

  // First collect existing leaves
  std::vector<BroadphaseObjectInfo<S>> leaf_nodes;
  for (const auto& kv : user_id_map_) {
    // Obtain the leaf node
    const auto leaf_node_index = kv.second;
    const auto& leaf_node = node_allocator_.Get(leaf_node_index);

    // Check the status of the leaf node
    assert(leaf_node.status.IsLeaf());
    if (leaf_node.status.IsRemoved()) {
      continue;
    }

    leaf_nodes.emplace_back(
        BroadphaseObjectInfo<S>{leaf_node.bv, leaf_node.user_id});
  }

  // Build the tree
  auto new_root_index = BuildTreeExternal(leaf_nodes.data(), leaf_nodes.size(),
                                          node_allocator_, &state.user_id_map);
  state.root_index = new_root_index;
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::ApplyUpdateStructure(
    TreeUpdateState&& state) {
  DestructTreeExternal(root_node_, node_allocator_);
  root_node_ = state.root_index;
  user_id_map_ = std::move(state.user_id_map);
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::AbortUpdateStructure(
    TreeUpdateState& state) {
  DestructTreeExternal(state.root_index, node_allocator_);
  state.user_id_map.clear();
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::PrepareAddNewObjects(
    BroadphaseObjectInfo<S>* new_objects, std::uint32_t n_new_objects,
    TreeUpdateState& state) {
  if (n_new_objects == 0U) {
    state.root_index = kInvalidAllocatorIndex;
    state.user_id_map.clear();
    return;
  }

  // Build the tree
  auto new_root_index = BuildTreeExternal(new_objects, n_new_objects,
                                          node_allocator_, &state.user_id_map);
  state.root_index = new_root_index;
}

template <typename S, template <typename Object> class ObjectAllocator>
bool BinaryAABB_Tree<S, ObjectAllocator>::ApplyAddNewObjects(
    TreeUpdateState& state) {
  // Special case 1: no new objects
  if (state.root_index == kInvalidAllocatorIndex) {
    return true;
  }

  // Special case 2: no existing objects
  if (root_node_ == kInvalidAllocatorIndex) {
    root_node_ = state.root_index;
    user_id_map_ = state.user_id_map;
    return true;
  }

  // General case
  const auto new_root_index = node_allocator_.AllocateObject();
  if (new_root_index == kInvalidAllocatorIndex) {
    return false;
  }

  // Obtain existing roots
  auto& old_root = node_allocator_.Get(root_node_);
  auto& inserted_root = node_allocator_.Get(state.root_index);

  // Assign meta
  auto& new_root = node_allocator_.Get(new_root_index);
  new_root.bv = old_root.bv + inserted_root.bv;
  new_root.children[0] = root_node_;
  new_root.children[1] = state.root_index;
  new_root.parent = kInvalidAllocatorIndex;
  new_root.status.SetAsInner();

  // Update the parent of old roots
  old_root.parent = new_root_index;
  inserted_root.parent = new_root_index;

  // Update user_id_map
  user_id_map_.insert(state.user_id_map.begin(), state.user_id_map.end());

  // Done
  return true;
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::AbortAddNewObjects(
    TreeUpdateState& state) {
  AbortUpdateStructure(state);
}

template <typename S, template <typename Object> class ObjectAllocator>
bool BinaryAABB_Tree<S, ObjectAllocator>::RemoveObject(
    std::uint64_t object_user_id) {
  // Find the leaf index from user map
  auto iter = user_id_map_.find(object_user_id);
  if (iter == user_id_map_.end()) {
    return false;
  }

  // Obtain the node
  const auto node_index = iter->second;
  auto& node = node_allocator_.Get(node_index);

  // Update AABB as empty (thus no collision can found it), and status as
  // removed (Actual removed happens at next structure update).
  node.bv = AABB<S>();
  node.status.SetAsRemoved();
  return true;
}

template <typename S, template <typename Object> class ObjectAllocator>
bool BinaryAABB_Tree<S, ObjectAllocator>::UpdateObjectAABB(
    std::uint64_t object_user_id, const AABB<S>& new_AABB) {
  // Find the leaf index from user map
  auto iter = user_id_map_.find(object_user_id);
  if (iter == user_id_map_.end()) {
    return false;
  }

  // Obtain the node
  auto node_index = iter->second;
  while (node_index != kInvalidAllocatorIndex) {
    // Copy the old bv
    auto& node = node_allocator_.Get(node_index);
    if (node.bv.contain(new_AABB)) {
      break;
    }

    // Assign the new bv and update to parent
    node.bv += new_AABB;
    node_index = node.parent;
  }

  // Done
  return true;
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::VisitTree(
    const VisitorFn& visitor_fn) const {
  // Special case of empty tree
  if (root_node_ == kInvalidAllocatorIndex) {
    return;
  }

  // Non-empty case
  std::stack<AllocatorIndex> task_stack;
  task_stack.push(root_node_);
  while (!task_stack.empty()) {
    // Obtain the node
    const auto node_index = task_stack.top();
    task_stack.pop();
    assert(node_index != kInvalidAllocatorIndex);
    const auto& node = node_allocator_.Get(node_index);

    // Visit it
    const bool is_leaf = node.status.IsLeaf();
    bool node_finish = false;
    bool overall_finish = false;
    visitor_fn(node.bv, is_leaf, node.user_id, node_finish, overall_finish);

    // Processing output
    if (overall_finish) {
      return;
    }
    if (node_finish || is_leaf) {
      continue;
    }

    // Handle its children
    assert(node.status.IsInner());
    assert(node.children[1] != kInvalidAllocatorIndex);
    assert(node.children[0] != kInvalidAllocatorIndex);
    task_stack.push(node.children[1]);
    task_stack.push(node.children[0]);
  }
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::SingleObjectCollision(
    const BroadphaseObjectInfo<S>& object, const CollisionFn& collision_fn,
    void* collision_fn_data) const {
  // Special case of empty tree
  if (root_node_ == kInvalidAllocatorIndex) {
    return;
  }

  // Non-empty case
  const AABB<S>& object_aabb = object.bv;
  std::stack<AllocatorIndex> task_stack;
  task_stack.push(root_node_);
  while (!task_stack.empty()) {
    // Obtain the node
    const auto node_index = task_stack.top();
    task_stack.pop();
    assert(node_index != kInvalidAllocatorIndex);
    const Node node = node_allocator_.Get(node_index);

    // Check bv
    const AABB<S>& node_bv = node.bv;
    if (!object_aabb.overlap(node_bv)) {
      continue;
    }

    // A valid instance
    if (node.status.IsLeaf()) {
      const bool overall_done =
          collision_fn(node.user_id, object.user_id, collision_fn_data);
      if (overall_done) {
        return;
      }
    } else {
      assert(node.status.IsInner());
      assert(node.children[1] != kInvalidAllocatorIndex);
      assert(node.children[0] != kInvalidAllocatorIndex);
      task_stack.push(node.children[1]);
      task_stack.push(node.children[0]);
    }
  }
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::TreeCollision(
    const BinaryAABB_Tree<S, ObjectAllocator>& tree2,
    const CollisionFn& collision_fn, void* collision_fn_data) const {
  // Special case of empty tree
  if (root_node_ == kInvalidAllocatorIndex ||
      tree2.root_node_ == kInvalidAllocatorIndex) {
    return;
  }

  // Non-empty case
  std::vector<std::pair<AllocatorIndex, AllocatorIndex>> task_stack;
  task_stack.emplace_back(root_node_, tree2.root_node_);
  while (!task_stack.empty()) {
    // Obtain the index
    const auto this_task = task_stack.back();
    task_stack.pop_back();
    const auto node1_index = this_task.first;
    const auto node2_index = this_task.second;
    assert(node1_index != kInvalidAllocatorIndex);
    assert(node2_index != kInvalidAllocatorIndex);

    // Obtain the node
    const Node node1 = node_allocator_.Get(node1_index);
    const Node node2 = tree2.node_allocator_.Get(node2_index);
    if (!node1.bv.overlap(node2.bv)) {
      continue;
    }

    // Both are leaf
    const auto is_leaf_1 = node1.status.IsLeaf();
    const auto is_leaf_2 = node2.status.IsLeaf();
    if (is_leaf_1 && is_leaf_2) {
      const bool overall_done =
          collision_fn(node1.user_id, node2.user_id, collision_fn_data);
      if (overall_done) {
        return;
      }

      // Done with both leaf
      continue;
    }

    // At least one is not leaf
    if (is_leaf_2 || ((!is_leaf_1) && (node1.bv.size() > node2.bv.size()))) {
      assert(node1.status.IsInner());
      assert(node1.children[1] != kInvalidAllocatorIndex);
      assert(node1.children[0] != kInvalidAllocatorIndex);
      task_stack.emplace_back(node1.children[1], node2_index);
      task_stack.emplace_back(node1.children[0], node2_index);
    } else {
      assert(node2.status.IsInner());
      assert(node2.children[1] != kInvalidAllocatorIndex);
      assert(node2.children[0] != kInvalidAllocatorIndex);
      task_stack.emplace_back(node1_index, node2.children[1]);
      task_stack.emplace_back(node1_index, node2.children[0]);
    }
  }
}

template <typename S, template <typename Object> class ObjectAllocator>
void BinaryAABB_Tree<S, ObjectAllocator>::SelfCollision(
    const CollisionFn& collision_fn, void* collision_fn_data) const {
  if (root_node_ == kInvalidAllocatorIndex) {
    return;
  }

  // Processing loop
  std::stack<std::pair<AllocatorIndex, AllocatorIndex>> task_stack;
  task_stack.push(std::make_pair(root_node_, root_node_));
  while (!task_stack.empty()) {
    const auto this_task = task_stack.top();
    task_stack.pop();
    const auto node1_index = this_task.first;
    const auto node2_index = this_task.second;
    if (node1_index == node2_index) {
      // Push the (left, right) child pair into the node
      const auto& node = node_allocator_.Get(node1_index);
      if (node.status.IsInner()) {
        task_stack.push(std::make_pair(node.children[0], node.children[0]));
        task_stack.push(std::make_pair(node.children[1], node.children[1]));
        task_stack.push(std::make_pair(node.children[0], node.children[1]));
      }

      // Done
      continue;
    }

    // Case of different node
    assert(node1_index != node2_index);
    const Node node1 = node_allocator_.Get(node1_index);
    const Node node2 = node_allocator_.Get(node2_index);
    const AABB<S>& node1_bv = node1.bv;
    const AABB<S>& node2_bv = node2.bv;
    if (!node1_bv.overlap(node2_bv)) {
      continue;
    }

    // Both are leaf
    const bool is_node1_leaf = node1.status.IsLeaf();
    const bool is_node2_leaf = node2.status.IsLeaf();
    if (is_node1_leaf && is_node2_leaf) {
      const bool overall_done =
          collision_fn(node1.user_id, node2.user_id, collision_fn_data);
      if (overall_done) {
        return;
      }

      // Done with both leaf
      continue;
    }

    // At least one is not leaf
    const bool continue_on_node1 =
        is_node2_leaf ||
        ((!is_node1_leaf) && (node1_bv.size() > node2_bv.size()));
    if (continue_on_node1) {
      assert(node1.status.IsInner());
      assert(node1.children[1] != kInvalidAllocatorIndex);
      assert(node1.children[0] != kInvalidAllocatorIndex);
      task_stack.push(std::make_pair(node1.children[1], node2_index));
      task_stack.push(std::make_pair(node1.children[0], node2_index));
    } else {
      assert(node2.status.IsInner());
      assert(node2.children[1] != kInvalidAllocatorIndex);
      assert(node2.children[0] != kInvalidAllocatorIndex);
      task_stack.push(std::make_pair(node1_index, node2.children[1]));
      task_stack.push(std::make_pair(node1_index, node2.children[0]));
    }
  }
}

template <typename S, template <typename Object> class ObjectAllocator>
bool BinaryAABB_Tree<S, ObjectAllocator>::SanityCheck() const {
  if (root_node_ == kInvalidAllocatorIndex) {
    return true;
  }

  std::stack<AllocatorIndex> task_stack;
  task_stack.push(root_node_);
  while (!task_stack.empty()) {
    const auto node_index = task_stack.top();
    task_stack.pop();
    const Node& node = node_allocator_.Get(node_index);

    // Connection
    if (node.parent != kInvalidAllocatorIndex) {
      const Node& parent_node = node_allocator_.Get(node.parent);
      if (parent_node.status.IsRemoved() || parent_node.status.IsLeaf()) {
        return false;
      }
      if (!parent_node.bv.contain(node.bv)) {
        return false;
      }
      if ((parent_node.children[0] != node_index) &&
          (parent_node.children[1] != node_index)) {
        return false;
      }
    }

    // Push the children
    if (node.status.IsInner()) {
      task_stack.push(node.children[1]);
      task_stack.push(node.children[0]);
    }
  }

  return true;
}

}  // namespace detail
}  // namespace fcl