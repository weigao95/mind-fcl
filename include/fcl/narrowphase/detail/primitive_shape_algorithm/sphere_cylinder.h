/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2018, Toyota Research Institute.
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

/** @author Sean Curtis <sean@tri.global> (2018) */

#ifndef FCL_NARROWPHASE_DETAIL_SPHERECYLINDER_H
#define FCL_NARROWPHASE_DETAIL_SPHERECYLINDER_H

#include "fcl/geometry/shape/cylinder.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/narrowphase/contact_point.h"

namespace fcl {

namespace detail {

/** @name       Custom cylinder-sphere proximity algorithms

 These functions provide custom algorithms for analyzing the relationship
 between a sphere and a cylinder.

 These functions make use of the
 [Drake monogram
 notation](http://drake.mit.edu/doxygen_cxx/group__multibody__notation__basics.html)
 to describe quantities (particularly the poses of shapes).

 Both shapes must be posed in a common frame (notated as F). This common frame
 is typically the world frame W. Regardless, if the optional output data is
 returned, it will be reported in this common frame F.

 The functions can be executed in one of two ways: to perform a strict boolean
 query (is colliding/is separated) or to get data which characterizes the
 colliding or separating state (e.g., penetration depth vs separating distance).
 The difference in usage is defined by whether the optional out parameters
 are null or non-null. In the documentation, these are labeled "(optional)".

 For these functions, if the sphere and cylinder are detected to be *touching*
 this is considered a collision. As such, a collision query would report true
 and a separating distance query would report false.
 */

// NOTE: the choice to consider touching contact as collision is predicated on
// the implementation in sphere-sphere contact.

//@{

/** Detect collision between the sphere and cylinder. If colliding, return
 characterization of collision in the provided vector.

 The reported depth is guaranteed to be continuous with respect to the relative
 pose. In contrast, the normal and contact position are only guaranteed to be
 similarly continuous while the sphere center lies *outside* the cylinder.
 However, if the sphere center lies *inside* the cylinder, there are regions of
 discontinuity in both normal and contact position. This is due to the fact that
 there is not necessarily a *unique* characterization of penetration depth
 (e.g., the sphere center may be equidistant to both a cap face as well as the
 barrel). This ambiguity is resolved through an arbitrary prioritization scheme.
 If the sphere center is equidistant to both an end face and the barrel, the end
 face is used for normal and contact position computation. If the sphere center
 is closer to the barrel, but there is no unique solution (i.e., the sphere
 center lies on the center axis), the sphere is arbitrarily considered to be
 penetrating from the direction of the cylinder's +x direction (yielding a
 contact normal in the cylinder's -x direction.)

 @param sphere         The sphere geometry.
 @param X_FS           The pose of the sphere S in the common frame F.
 @param cylinder       The cylinder geometry.
 @param X_FC           The pose of the cylinder C in the common frame F.
 @param contacts[out]  (optional) If the shapes collide, the contact point data
                       will be appended to the end of this vector.
 @return True if the objects are colliding (including touching).
 @tparam S The scalar parameter (must be a valid Eigen scalar).  */
template <typename S>
bool sphereCylinderIntersect(const Sphere<S>& sphere, const Transform3<S>& X_FS,
                             const Cylinder<S>& cylinder,
                             const Transform3<S>& X_FC,
                             CollisionPenetrationContactData<S>* contacts);

/** Evaluate the minimum separating distance between a sphere and cylinder. If
 separated, the nearest points on each shape will be returned in frame F.
 @param sphere         The sphere geometry.
 @param X_FS           The pose of the sphere S in the common frame F.
 @param cylinder       The cylinder geometry.
 @param X_FC           The pose of the cylinder C in the common frame F.
 @param distance[out]  (optional) The separating distance between the cylinder
                       and sphere. Set to -1 if the shapes are penetrating.
 @param p_FSc[out]     (optional) The closest point on the *sphere* to the
                       cylinder measured and expressed in frame F.
 @param p_FCs[out]     (optional) The closest point on the *cylinder* to the
                       sphere measured and expressed in frame F.
 @return True if the objects are separated.
 @tparam S The scalar parameter (must be a valid Eigen scalar).  */
template <typename S>
bool sphereCylinderDistance(const Sphere<S>& sphere, const Transform3<S>& X_FS,
                            const Cylinder<S>& cylinder,
                            const Transform3<S>& X_FC, S* distance,
                            Vector3<S>* p_FSc, Vector3<S>* p_FCs);

//@}

}  // namespace detail
}  // namespace fcl

#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_cylinder-inl.h"

#endif  // FCL_NARROWPHASE_DETAIL_SPHERECYLINDER_H
