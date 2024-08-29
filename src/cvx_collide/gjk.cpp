#include "fcl/cvx_collide/gjk.h"

namespace fcl {
namespace cvx_collide {

template struct GJKSimplex<float>;
template struct GJKSimplex<double>;

template class GJK<float>;
template class GJK<double>;

}  // namespace cvx_collide
}  // namespace fcl