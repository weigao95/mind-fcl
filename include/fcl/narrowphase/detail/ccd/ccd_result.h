//
// Created by Wei Gao on 2024/6/13.
//

#pragma once

#include "fcl/narrowphase/detail/ccd/ccd_contact.h"

namespace fcl {

template <typename S>
struct ContinuousCollisionResult {
 public:
  using UserContactProcessFunctor = std::function<void(
      const ContinuousCollisionContact<S>& c, bool& output_keep_this_contact,
      bool& output_can_we_terminate)>;
  explicit ContinuousCollisionResult();
  explicit ContinuousCollisionResult(UserContactProcessFunctor user_fn);
  ~ContinuousCollisionResult() = default;

  /// The interface for contact
  void AddContact(const ContinuousCollisionContact<S>& c);
  void ClearContact();
  std::size_t num_contacts() const;
  const std::vector<ContinuousCollisionContact<S>>& raw_contacts() const;

  /// Determine whether the collision checking can stop
  bool TerminationConditionSatisfied(std::size_t n_max_contacts) const;

 private:
  // Raw contact data
  std::vector<ContinuousCollisionContact<S>> contacts_;

  // Optional user specified functor
  UserContactProcessFunctor user_process_functor_;
  bool user_stop_{false};
};

}  // namespace fcl

#include "fcl/narrowphase/detail/ccd/ccd_result-inl.h"