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

#include "test_fcl_utility.h"

#include <cstddef>
#include <cstdio>
#include <unordered_map>

#include "fcl/narrowphase/collision.h"

namespace fcl {

namespace test {

//==============================================================================
Timer::Timer() {
#ifdef _WIN32
  QueryPerformanceFrequency(&frequency);
  startCount.QuadPart = 0;
  endCount.QuadPart = 0;
#else
  startCount.tv_sec = startCount.tv_usec = 0;
  endCount.tv_sec = endCount.tv_usec = 0;
#endif

  stopped = 0;
  startTimeInMicroSec = 0;
  endTimeInMicroSec = 0;
}

//==============================================================================
Timer::~Timer() {
  // Do nothing
}

//==============================================================================
void Timer::start() {
  stopped = 0;  // reset stop flag
#ifdef _WIN32
  QueryPerformanceCounter(&startCount);
#else
  gettimeofday(&startCount, nullptr);
#endif
}

//==============================================================================
void Timer::stop() {
  stopped = 1;  // set timer stopped flag

#ifdef _WIN32
  QueryPerformanceCounter(&endCount);
#else
  gettimeofday(&endCount, nullptr);
#endif
}

double Timer::getElapsedTimeInMicroSec() {
#ifdef _WIN32
  if (!stopped) QueryPerformanceCounter(&endCount);

  startTimeInMicroSec = startCount.QuadPart * (1000000.0 / frequency.QuadPart);
  endTimeInMicroSec = endCount.QuadPart * (1000000.0 / frequency.QuadPart);
#else
  if (!stopped) gettimeofday(&endCount, nullptr);

  startTimeInMicroSec = (startCount.tv_sec * 1000000.0) + startCount.tv_usec;
  endTimeInMicroSec = (endCount.tv_sec * 1000000.0) + endCount.tv_usec;
#endif

  return endTimeInMicroSec - startTimeInMicroSec;
}

//==============================================================================
double Timer::getElapsedTimeInMilliSec() {
  return this->getElapsedTimeInMicroSec() * 0.001;
}

//==============================================================================
double Timer::getElapsedTimeInSec() {
  return this->getElapsedTimeInMicroSec() * 0.000001;
}

//==============================================================================
double Timer::getElapsedTime() { return this->getElapsedTimeInMilliSec(); }

//==============================================================================
std::string getNodeTypeName(NODE_TYPE node_type) {
  if (node_type == BV_UNKNOWN)
    return std::string("BV_UNKNOWN");
  else if (node_type == BV_AABB)
    return std::string("BV_AABB");
  else if (node_type == BV_OBB)
    return std::string("BV_OBB");
  else if (node_type == BV_RSS)
    return std::string("BV_RSS");
  else if (node_type == BV_kIOS)
    return std::string("BV_kIOS");
  else if (node_type == BV_OBBRSS)
    return std::string("BV_OBBRSS");
  else if (node_type == BV_KDOP16)
    return std::string("BV_KDOP16");
  else if (node_type == BV_KDOP18)
    return std::string("BV_KDOP18");
  else if (node_type == BV_KDOP24)
    return std::string("BV_KDOP24");
  else if (node_type == GEOM_BOX)
    return std::string("GEOM_BOX");
  else if (node_type == GEOM_SPHERE)
    return std::string("GEOM_SPHERE");
  else if (node_type == GEOM_ELLIPSOID)
    return std::string("GEOM_ELLIPSOID");
  else if (node_type == GEOM_CAPSULE)
    return std::string("GEOM_CAPSULE");
  else if (node_type == GEOM_CONE)
    return std::string("GEOM_CONE");
  else if (node_type == GEOM_CYLINDER)
    return std::string("GEOM_CYLINDER");
  else if (node_type == GEOM_CONVEX)
    return std::string("GEOM_CONVEX");
  else if (node_type == GEOM_PLANE)
    return std::string("GEOM_PLANE");
  else if (node_type == GEOM_HALFSPACE)
    return std::string("GEOM_HALFSPACE");
  else if (node_type == GEOM_TRIANGLE)
    return std::string("GEOM_TRIANGLE");
  else if (node_type == GEOM_OCTREE2)
    return std::string("GEOM_OCTREE2");
  else
    return std::string("invalid");
}

//==============================================================================
bool getTestModelDirectory(std::string& test_model_dir) {
#ifdef PRIVATE_TEST_MODEL_DIR
  auto filepath = std::string(PRIVATE_TEST_MODEL_DIR);
  test_model_dir = filepath;
  return true;
#else
  FCL_UNUSED(test_model_dir);
  return false;
#endif
}

std::shared_ptr<const octomap::Pointcloud> readPointCloudXYZ(
    const std::string& pcl_xyz_path) {
  // Read the cloud
  std::ifstream read_stream(pcl_xyz_path, std::ios_base::in);
  std::string line;
  auto point_cloud = std::make_shared<octomap::Pointcloud>();

  while (std::getline(read_stream, line)) {
    // Process str
    // std::cout << line << std::endl;
    std::stringstream line_ss(line);
    float x, y, z;
    line_ss >> x;
    line_ss >> y;
    line_ss >> z;

    // Make the cloud
    octomap::point3d this_point(x, y, z);
    point_cloud->push_back(this_point);
  }

  // Done
  read_stream.close();
  return point_cloud;
}

std::shared_ptr<Octree2CollisionGeometry<float>> buildOctree2(
    const std::string& pcl_xyz_path, float voxel_resolution,
    std::uint16_t bottom_half_shape) {
  // Read the cloud
  auto point_cloud = readPointCloudXYZ(pcl_xyz_path);
  if (point_cloud == nullptr) return nullptr;

  // Point fn
  auto point_fn = [&point_cloud](int index, float& x, float& y,
                                 float& z) -> void {
    const auto point = point_cloud->getPoint(index);
    x = point.x();
    y = point.y();
    z = point.z();
  };

  auto tree = std::make_shared<octree2::Octree<float>>(voxel_resolution,
                                                       bottom_half_shape);
  tree->rebuildTree(point_fn, point_cloud->size());
  auto tree_geom = std::make_shared<Octree2CollisionGeometry<float>>(tree);
  tree_geom->computeLocalAABB();
  return tree_geom;
}

}  // namespace test
}  // namespace fcl

/// Load stl mesh
using FlatPoint = std::array<float, 3>;
template <>
struct std::equal_to<FlatPoint> {
  bool operator()(const FlatPoint& left, const FlatPoint& right) const {
    return left[0] == right[0] && left[1] == right[1] && left[2] == right[2];
  }
};

struct Hasher {
  std::size_t operator()(const FlatPoint& point) const {
    return std::hash<float>()(point[0]) * std::hash<float>()(point[1]) ^
           std::hash<float>()(point[2]);
  }
};

FlatPoint readQVec3DFromBinaryStream(std::ifstream& ifs) {
  FlatPoint vec;
  for (int i = 0; i < 3; ++i) {
    float f;
    ifs.read(reinterpret_cast<char*>(&f), sizeof(float));
    vec[i] = f;
  }
  return vec;
}

void fcl::test::loadSTLMesh(const std::string& mesh_path,
                            std::vector<fcl::Vector3f>& points,
                            std::vector<fcl::MeshSimplex>& triangles) {
  points.clear();
  triangles.clear();
  std::ifstream ifs(mesh_path, std::ios::in | std::ios::binary);

  constexpr std::streamsize header_size = 80;
  constexpr std::streamsize facets_size = 50;

  // Read header
  char header[header_size];
  ifs.read(header, header_size);
  if (ifs.gcount() != header_size) {
    throw std::runtime_error("Header read failure");
  }

  // Read facets number
  unsigned facets_num = 0;
  ifs.read(reinterpret_cast<char*>(&facets_num), 4);

  // Check file size
  const auto current = ifs.tellg();
  ifs.seekg(0, std::ios::end);
  if (facets_num * facets_size != ifs.tellg() - current) {
    throw std::runtime_error("Header read failure");
  }
  ifs.seekg(current, std::ios::beg);

  std::unordered_map<FlatPoint, unsigned, Hasher> pointIndices;
  unsigned facetCount = facets_num;
  unsigned vertexCount = 0;
  constexpr std::streamsize byteLengthForVertex = 3 * sizeof(float);
  constexpr std::streamsize byteLengthForAttr = 2;
  while (facetCount-- != 0) {
    char buf[byteLengthForVertex];

    // read normal
    ifs.read(buf, byteLengthForVertex);

    // read vertices
    std::array<unsigned, 3> flat_triangle;
    for (unsigned ots = 0; ots < 3; ++ots) {
      const FlatPoint newPoint = readQVec3DFromBinaryStream(ifs);
      auto item = pointIndices.emplace(newPoint, vertexCount);
      if (item.second) {
        points.push_back(fcl::Vector3f(newPoint[0], newPoint[1], newPoint[2]));
        ++vertexCount;
      }
      flat_triangle.at(ots) = item.first->second;
    }

    fcl::MeshSimplex fcl_triangle;
    fcl_triangle.set(flat_triangle[0], flat_triangle[1], flat_triangle[2]);
    triangles.push_back(fcl_triangle);

    // read attributes
    ifs.read(buf, byteLengthForAttr);
  }
}

struct BinvoxLoadingPara {
  BinvoxLoadingPara() = default;
  int version = 1;           // binvox file format version (should be 1)
  int depth, height, width;  // dimensions of the voxel grid
  int size;                  // number of grid cells (height * width * depth)
  double tx, ty, tz;         // Translation
  double scale;              // Scaling factor
  double resolution;         // resolution for the octree
};

std::string readBinvoxHeader(std::ifstream& ifs, BinvoxLoadingPara& para) {
  using namespace std;
  string line;
  ifs >> line;
  if (line.compare("#binvox") != 0) {
    return "First line is not '#binvox'";
  }
  ifs >> para.version;
  para.depth = -1;
  bool done = false;
  while (ifs.good() && !done) {
    ifs >> line;
    if (line.compare("data") == 0)
      done = true;
    else if (line.compare("dim") == 0) {
      ifs >> para.depth >> para.height >> para.width;
      para.size = para.depth * para.height * para.width;
    } else if (line.compare("translate") == 0) {
      ifs >> para.tx >> para.ty >> para.tz;
    } else if (line.compare("scale") == 0) {
      ifs >> para.scale;
    } else {
      char c;
      do {  // skip until end of line
        c = ifs.get();
      } while (ifs.good() && (c != '\n'));
    }
  }
  if (!done) {
    return "Error reading binvox header.";
  }
  if (para.depth == -1) {
    return "Error missing dimensions in header while reading binvox.";
  }
  const int maxSide = std::max(std::max(para.width, para.height), para.depth);
  para.resolution = double(para.scale) / double(maxSide);
  return {};
}

std::shared_ptr<fcl::octree2::Octree<float>> fcl::test::loadBinvoxAsOctree2(
    const std::string& filepath, double scale,
    std::uint16_t bottom_half_shape) {
  using namespace std;
  ifstream ifs(filepath, ios::in | ios::binary);
  if (!ifs.good()) {
    return nullptr;
  }
  BinvoxLoadingPara para;
  auto errMsg = readBinvoxHeader(ifs, para);
  if (!errMsg.empty()) {
    std::cout << errMsg;
    return nullptr;
  }
  para.tx *= scale;
  para.ty *= scale;
  para.tz *= scale;
  para.resolution *= scale;
  unsigned char value;
  unsigned char count;
  int index = 0;
  int end_index = 0;
  // unsigned nr_voxels = 0;
  ifs.unsetf(ios::skipws);  // need to read every byte now (!)
  ifs >> value;             // read the linefeed char
  fcl::octomap::Pointcloud loaded_point;
  while ((end_index < para.size) && ifs.good()) {
    ifs >> value >> count;
    if (ifs.good()) {
      end_index = index + count;
      if (end_index > para.size) return nullptr;
      for (int j = index; j < end_index; j++) {
        // voxel index --> voxel coordinates
        const int y = j % para.width;
        const int z = (j / para.width) % para.height;
        const int x = j / (para.width * para.height);

        // voxel coordinates --> world coordinates
        double point_x = (x * para.resolution + para.tx);
        double point_y = (y * para.resolution + para.ty);
        double point_z = (z * para.resolution + para.tz);
        if (value == 1) {
          loaded_point.push_back(point_x, point_y, point_z);
        }
      }

      // if (value > 0) nr_voxels += count;
      index = end_index;
    }  // if file still ok
  }    // while
  ifs.close();

  auto point_fn = [&loaded_point](int index, float& x, float& y,
                                  float& z) -> void {
    const auto point = loaded_point[index];
    x = point.x();
    y = point.y();
    z = point.z();
  };
  auto tree = std::make_shared<fcl::octree2::Octree<float>>(para.resolution,
                                                            bottom_half_shape);
  tree->rebuildTree(point_fn, loaded_point.size());
  return tree;
}