#include <vector>
#include <chrono>
#include "fcl/geometry/heightmap/flat_heightmap.h"

namespace fcl {
namespace heightmap {

std::shared_ptr<const PointCloud> generateRandomPointCloud(
    std::size_t numPoints) {
  auto pointCloud = std::make_shared<PointCloud>();
  pointCloud->reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Eigen::Vector3f point = Eigen::Vector3f::Random();
    pointCloud->push_back(point[0], point[1], point[2] + 1);
  }
  return pointCloud;
}

template <typename S>
void rebuildMapBenchmark() {
  for (std::size_t size : {10000, 100000, 1000000}) {
    const std::shared_ptr<const fcl::octomap::Pointcloud> points =
        generateRandomPointCloud(size);
    for (int half_map_shape = 64; half_map_shape <= 512; half_map_shape *= 2) {
      for (float resolution = 0.001; resolution < 0.01; resolution *= 2) {
        auto start = std::chrono::high_resolution_clock::now();
        auto heightMap = std::make_shared<FlatHeightMap<S>>(
            resolution, half_map_shape);
        heightMap->resetHeights();
        heightMap->updateHeightsByPointCloud3D(*points);
        auto end = std::chrono::high_resolution_clock::now();
        auto ms_time =
            std::chrono::duration_cast<std::chrono::milliseconds>((end - start))
                .count();
        std::cout << "rebuildMap Time in ms: " << ms_time
                  << " size of pointcloud: " << size
                  << " half_map_shape: " << half_map_shape
                  << " resolution: " << resolution << std::endl;
      }
    }
  }
}

}  // namespace heightmap
}  // namespace fcl


//==============================================================================
int main() {
  std::cout << "Benchmark with float" << std::endl;
  fcl::heightmap::rebuildMapBenchmark<float>();
  std::cout << "Benchmark with double" << std::endl;
  fcl::heightmap::rebuildMapBenchmark<double>();
}
