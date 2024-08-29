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

#ifndef FCL_COLLISION_OBJECT_INL_H
#define FCL_COLLISION_OBJECT_INL_H

#include "fcl/narrowphase/collision_object.h"

namespace fcl
{
//==============================================================================
template <typename S>
CollisionObject<S>::CollisionObject(
    const std::shared_ptr<CollisionGeometry<S>>& cgeom_)
  : cgeom_const(cgeom_), t(Transform3<S>::Identity())
{
  if (cgeom_ != nullptr)
  {
    cgeom_->computeLocalAABB();
    computeAABB();
  }
}

//==============================================================================
template <typename S>
CollisionObject<S>::CollisionObject(
    const std::shared_ptr<CollisionGeometry<S>>& cgeom_,
    const Transform3<S>& tf)
  : cgeom_const(cgeom_), t(tf)
{
  cgeom_->computeLocalAABB();
  computeAABB();
}

//==============================================================================
template <typename S>
CollisionObject<S>::CollisionObject(
    const std::shared_ptr<CollisionGeometry<S>>& cgeom_,
    const Matrix3<S>& R,
    const Vector3<S>& T)
  : cgeom_const(cgeom_), t(Transform3<S>::Identity())
{
  t.linear() = R;
  t.translation() = T;
  cgeom_->computeLocalAABB();
  computeAABB();
}

//==============================================================================
template <typename S>
CollisionObject<S>::~CollisionObject()
{
  // Do nothing
}

//==============================================================================
template <typename S>
CollisionObject<S>::CollisionObject(
    const std::shared_ptr<const CollisionGeometry<S>>& cgeom,
    const GeometryLocalAABBComputed&)
    : cgeom_const(cgeom), t(Transform3<S>::Identity())
{
  if (cgeom != nullptr)
  {
    // Should be invoked externally
    // cgeom_->computeLocalAABB();
    computeAABB();
  }
}

//==============================================================================
template <typename S>
CollisionObject<S>::CollisionObject(
    const std::shared_ptr<const CollisionGeometry<S>>& cgeom,
    const Transform3<S>& tf,
    const GeometryLocalAABBComputed&)
    : cgeom_const(cgeom), t(tf)
{
  // Should be invoked externally
  // cgeom_->computeLocalAABB();
  computeAABB();
}

//==============================================================================
template <typename S>
OBJECT_TYPE CollisionObject<S>::getObjectType() const
{
  return cgeom_const->getObjectType();
}

//==============================================================================
template <typename S>
NODE_TYPE CollisionObject<S>::getNodeType() const
{
  return cgeom_const->getNodeType();
}

//==============================================================================
template <typename S>
const AABB<S>& CollisionObject<S>::getAABB() const
{
  return aabb;
}

//==============================================================================
template <typename S>
void CollisionObject<S>::computeAABB()
{
  if(t.linear().isIdentity())
  {
    aabb = translate(cgeom_const->aabb_local, t.translation());
  }
  else
  {
    Vector3<S> center = t * cgeom_const->aabb_center;
    Vector3<S> delta = Vector3<S>::Constant(cgeom_const->aabb_radius);
    aabb.min_ = center - delta;
    aabb.max_ = center + delta;
  }
}

//==============================================================================
template <typename S>
void*CollisionObject<S>::getUserData() const
{
  return user_data;
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setUserData(void* data)
{
  user_data = data;
}

//==============================================================================
template <typename S>
const Vector3<S> CollisionObject<S>::getTranslation() const
{
  return t.translation();
}

//==============================================================================
template <typename S>
const Matrix3<S> CollisionObject<S>::getRotation() const
{
  return t.linear();
}

//==============================================================================
template <typename S>
const Quaternion<S> CollisionObject<S>::getQuatRotation() const
{
  return Quaternion<S>(t.linear());
}

//==============================================================================
template <typename S>
const Transform3<S>& CollisionObject<S>::getTransform() const
{
  return t;
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setRotation(const Matrix3<S>& R)
{
  t.linear() = R;
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setTranslation(const Vector3<S>& T)
{
  t.translation() = T;
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setQuatRotation(const Quaternion<S>& q)
{
  t.linear() = q.toRotationMatrix();
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setTransform(const Matrix3<S>& R, const Vector3<S>& T)
{
  setRotation(R);
  setTranslation(T);
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setTransform(const Quaternion<S>& q, const Vector3<S>& T)
{
  setQuatRotation(q);
  setTranslation(T);
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setTransform(const Transform3<S>& tf)
{
  t = tf;
}

//==============================================================================
template <typename S>
bool CollisionObject<S>::isIdentityTransform() const
{
  return t.matrix().isIdentity();
}

//==============================================================================
template <typename S>
void CollisionObject<S>::setIdentityTransform()
{
  t.setIdentity();
}

//==============================================================================
template <typename S>
const std::shared_ptr<const CollisionGeometry<S>>&
CollisionObject<S>::collisionGeometry() const
{
  return cgeom_const;
}

//==============================================================================
template <typename S>
S CollisionObject<S>::getCostDensity() const
{
  return cgeom_const->cost_density;
}

//==============================================================================
template <typename S>
bool CollisionObject<S>::isOccupied() const
{
  return cgeom_const->isOccupied();
}

//==============================================================================
template <typename S>
bool CollisionObject<S>::isFree() const
{
  return cgeom_const->isFree();
}

//==============================================================================
template <typename S>
bool CollisionObject<S>::isUncertain() const
{
  return cgeom_const->isUncertain();
}

} // namespace fcl

#endif
