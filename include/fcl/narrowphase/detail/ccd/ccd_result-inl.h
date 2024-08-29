#pragma once

namespace fcl {

template <typename S>
ContinuousCollisionResult<S>::ContinuousCollisionResult() : user_stop_(false) {}

template <typename S>
ContinuousCollisionResult<S>::ContinuousCollisionResult(
    UserContactProcessFunctor user_fn)
    : user_process_functor_(std::move(user_fn)), user_stop_(false) {}

template <typename S>
void ContinuousCollisionResult<S>::AddContact(
    const ContinuousCollisionContact<S>& c) {
  if (user_process_functor_) {
    bool suggest_stop = false;
    bool keep_this = true;
    user_process_functor_(c, keep_this, suggest_stop);
    user_stop_ |= suggest_stop;
    if (keep_this) {
      contacts_.push_back(c);
    }
  } else {
    contacts_.emplace_back(c);
  }
}

template <typename S>
void ContinuousCollisionResult<S>::ClearContact() {
  contacts_.clear();
}

template <typename S>
std::size_t ContinuousCollisionResult<S>::num_contacts() const {
  return contacts_.size();
}

template <typename S>
const std::vector<ContinuousCollisionContact<S>>&
ContinuousCollisionResult<S>::raw_contacts() const {
  return contacts_;
}

template <typename S>
bool ContinuousCollisionResult<S>::TerminationConditionSatisfied(
    std::size_t n_max_contacts) const {
  return user_stop_ || (contacts_.size() >= n_max_contacts);
}

}  // namespace fcl