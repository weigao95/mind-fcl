/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan */

#ifndef FCL_TRAVERSAL_RECURSE_INL_H
#define FCL_TRAVERSAL_RECURSE_INL_H

#include <queue>

#include "fcl/common/unused.h"
#include "fcl/narrowphase/detail/traversal/traversal_recurse.h"

namespace fcl {

namespace detail {

//==============================================================================
template <typename S>
void collisionRecurse(CollisionTraversalNodeBase<S>* node, int b1, int b2,
                      BVHFrontList* front_list) {
  const bool is_leaf_1 = node->isFirstNodeLeaf(b1);
  const bool is_leaf_2 = node->isSecondNodeLeaf(b2);

  if (is_leaf_1 && is_leaf_2) {
    updateFrontList(front_list, b1, b2);
    if (node->BVTesting(b1, b2)) return;
    node->leafTesting(b1, b2);
    return;
  }

  if (node->BVTesting(b1, b2)) {
    updateFrontList(front_list, b1, b2);
    return;
  }

  if (node->firstOverSecond(b1, b2)) {
    const int child_1 = node->getFirstLeftChild(b1);
    const int child_2 = node->getFirstRightChild(b1);
    collisionRecurse(node, child_1, b2, front_list);

    // early stop is disabled is front_list is used
    if (node->canStop() && !front_list) return;

    collisionRecurse(node, child_2, b2, front_list);
  } else {
    const int child_1 = node->getSecondLeftChild(b2);
    const int child_2 = node->getSecondRightChild(b2);
    collisionRecurse(node, b1, child_1, front_list);

    // early stop is disabled is front_list is used
    if (node->canStop() && !front_list) return;

    collisionRecurse(node, b1, child_2, front_list);
  }
}

//==============================================================================
/** @brief Bounding volume test structure */
template <typename S>
struct BVT {
  /** @brief distance between bvs */
  S d;

  /** @brief bv indices for a pair of bvs in two models */
  int b1, b2;
};

//==============================================================================
/** @brief Comparer between two BVT */
template <typename S>
struct BVT_Comparer {
  bool operator()(const BVT<S>& lhs, const BVT<S>& rhs) const {
    return lhs.d > rhs.d;
  }
};

//==============================================================================
template <typename S>
struct BVTQ {
  BVTQ() : qsize(2) {}

  bool empty() const { return pq.empty(); }

  size_t size() const { return pq.size(); }

  const BVT<S>& top() const { return pq.top(); }

  void push(const BVT<S>& x) { pq.push(x); }

  void pop() { pq.pop(); }

  bool full() const { return (pq.size() + 1 >= qsize); }

  std::priority_queue<BVT<S>, std::vector<BVT<S>>, BVT_Comparer<S>> pq;

  /** @brief Queue size */
  unsigned int qsize;
};

//==============================================================================
template <typename S>
void propagateBVHFrontListCollisionRecurse(CollisionTraversalNodeBase<S>* node,
                                           BVHFrontList* front_list) {
  BVHFrontList::iterator front_iter;
  BVHFrontList append;
  for (front_iter = front_list->begin(); front_iter != front_list->end();
       ++front_iter) {
    int b1 = front_iter->left;
    int b2 = front_iter->right;
    bool l1 = node->isFirstNodeLeaf(b1);
    bool l2 = node->isSecondNodeLeaf(b2);

    if (l1 & l2) {
      front_iter->valid = false;  // the front node is no longer valid, in
                                  // collideRecurse will add again.
      collisionRecurse(node, b1, b2, &append);
    } else {
      if (!node->BVTesting(b1, b2)) {
        front_iter->valid = false;

        if (node->firstOverSecond(b1, b2)) {
          int c1 = node->getFirstLeftChild(b1);
          int c2 = node->getFirstRightChild(b1);

          collisionRecurse(node, c1, b2, front_list);
          collisionRecurse(node, c2, b2, front_list);
        } else {
          int c1 = node->getSecondLeftChild(b2);
          int c2 = node->getSecondRightChild(b2);

          collisionRecurse(node, b1, c1, front_list);
          collisionRecurse(node, b1, c2, front_list);
        }
      }
    }
  }

  // clean the old front list (remove invalid node)
  for (front_iter = front_list->begin(); front_iter != front_list->end();) {
    if (!front_iter->valid)
      front_iter = front_list->erase(front_iter);
    else
      ++front_iter;
  }

  for (front_iter = append.begin(); front_iter != append.end(); ++front_iter) {
    front_list->push_back(*front_iter);
  }
}

}  // namespace detail
}  // namespace fcl

#endif
