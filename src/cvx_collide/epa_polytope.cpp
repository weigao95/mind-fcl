#include "fcl/cvx_collide/epa_polytope.h"

namespace fcl {
namespace cvx_collide {

template struct PolytopeVertex<float>;
template struct PolytopeVertex<double>;

template struct PolytopeEdge<float>;
template struct PolytopeEdge<double>;

template struct PolytopeFace<float>;
template struct PolytopeFace<double>;

template class Polytope<float>;
template class Polytope<double>;

}  // namespace cvx_collide
}  // namespace fcl