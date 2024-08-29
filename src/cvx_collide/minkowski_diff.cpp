#include "fcl/cvx_collide/minkowski_diff.h"

namespace fcl {
namespace cvx_collide {

template struct MinkowskiDiffVertex<float>;
template struct MinkowskiDiffVertex<double>;

template struct MinkowskiDiff<float>;
template struct MinkowskiDiff<double>;

}  // namespace cvx_collide
}  // namespace fcl