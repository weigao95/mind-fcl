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

#ifndef FCL_BROAD_PHASE_BRUTE_FORCE_INL_H
#define FCL_BROAD_PHASE_BRUTE_FORCE_INL_H

#include <iterator>

namespace fcl {

//==============================================================================
template <typename S>
NaiveCollisionManager<S>::NaiveCollisionManager()
{
  // Do nothing
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::registerObjects(const std::vector<CollisionObject<S>*>& other_objs)
{
  std::copy(other_objs.begin(), other_objs.end(), std::back_inserter(objs));
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::unregisterObject(CollisionObject<S>* obj)
{
  objs.remove(obj);
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::registerObject(CollisionObject<S>* obj)
{
  objs.push_back(obj);
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::setup()
{
  // Do nothing
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::update()
{
  // Do nothing
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::clear()
{
  objs.clear();
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::getObjects(std::vector<CollisionObject<S>*>& objs_) const
{
  objs_.resize(objs.size());
  std::copy(objs.begin(), objs.end(), objs_.begin());
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::collide(CollisionObject<S>* obj, void* cdata, CollisionCallBack<S> callback) const
{
  if(size() == 0) return;

  for(auto* obj2 : objs)
  {
    if(callback(obj, obj2, cdata))
      return;
  }
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::collide(void* cdata, CollisionCallBack<S> callback) const
{
  if(size() == 0) return;

  for(typename std::list<CollisionObject<S>*>::const_iterator it1 = objs.begin(), end = objs.end();
      it1 != end; ++it1)
  {
    typename std::list<CollisionObject<S>*>::const_iterator it2 = it1; it2++;
    for(; it2 != end; ++it2)
    {
      if((*it1)->getAABB().overlap((*it2)->getAABB()))
      {
        if(callback(*it1, *it2, cdata))
          return;
      }
    }
  }
}

//==============================================================================
template <typename S>
void NaiveCollisionManager<S>::collide(NaiveCollisionManager<S>* other_manager, void* cdata, CollisionCallBack<S> callback) const
{
  if((size() == 0) || (other_manager->size() == 0)) return;

  if(this == other_manager)
  {
    collide(cdata, callback);
    return;
  }

  for(auto* obj1 : objs)
  {
    for(auto* obj2 : other_manager->objs)
    {
      if(obj1->getAABB().overlap(obj2->getAABB()))
      {
        if(callback(obj1, obj2, cdata))
          return;
      }
    }
  }
}

//==============================================================================
template <typename S>
bool NaiveCollisionManager<S>::empty() const
{
  return objs.empty();
}

//==============================================================================
template <typename S>
size_t NaiveCollisionManager<S>::size() const
{
  return objs.size();
}

} // namespace fcl

#endif
