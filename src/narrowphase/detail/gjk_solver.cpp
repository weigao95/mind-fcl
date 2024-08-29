#include "fcl/narrowphase/detail/gjk_solver.h"

namespace fcl {
namespace detail {

template struct GJKSolver<float>;
template struct GJKSolver<double>;

}  // namespace detail
}  // namespace fcl