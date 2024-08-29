#===============================================================================
# Google Test settings
#===============================================================================
include_directories (BEFORE
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Build gtest
add_library(gtest STATIC gtest/src/gtest-all.cc)
add_library(gtest_main STATIC gtest/src/gtest_main.cc)
target_link_libraries(gtest_main gtest)
if(NOT WIN32)
  target_link_libraries(gtest pthread)
endif()

#execute_process(COMMAND cmake -E remove_directory ${CMAKE_BINARY_DIR}/test_results)
#execute_process(COMMAND cmake -E make_directory ${CMAKE_BINARY_DIR}/test_results)
include_directories(${GTEST_INCLUDE_DIRS})

add_library(test_fcl_utility test_fcl_utility.cpp profiler.cpp)
target_link_libraries(test_fcl_utility PUBLIC fcl)

# The internal test model
if (PRIVATE_TEST_MODEL_DIR)
    add_definitions(-DPRIVATE_TEST_MODEL_DIR="${PRIVATE_TEST_MODEL_DIR}")
endif ()

# test file list
set(tests
    test_fcl_box_box.cpp
    test_fcl_bvh_models.cpp
    test_fcl_bvh_models_tetrahedron.cpp
    test_fcl_tetrahedron_bvh_models_collision.cpp
    test_fcl_capsule_capsule.cpp
    test_fcl_cylinder_half_space.cpp
    test_fcl_collision.cpp
    test_fcl_collision_penetration.cpp
    test_fcl_constant_eps.cpp
    test_fcl_frontlist.cpp
    test_fcl_general.cpp
    test_fcl_generate_bvh_model_deferred_finalize.cpp
    test_fcl_geometric_shapes.cpp
    test_fcl_math.cpp
    test_fcl_profiler.cpp
    test_fcl_shape_mesh_consistency.cpp
    test_fcl_simple.cpp
    test_fcl_sphere_box.cpp
    test_fcl_sphere_capsule.cpp
    test_fcl_sphere_cylinder.cpp
    test_fcl_convex_primitive.cpp
    test_octree_from_triangles.cpp
    test_fixed_rotation_obb.cpp
    # gjk, epa and mpr
    cvx_collide/test_epa2.cpp
    cvx_collide/test_gjk2.cpp
    cvx_collide/test_gjk2_distance.cpp
    cvx_collide/test_epa2_with_gjk2.cpp
    cvx_collide/test_mpr_primitive.cpp
    cvx_collide/test_mpr_triangle.cpp
    cvx_collide/test_mpr_convex.cpp
    cvx_collide/test_mpr_penetration.cpp
    cvx_collide/test_mpr_local_refinement.cpp
    # broadphase
    broadphase/test_binary_AABB_tree_allocator.cpp
    broadphase/test_binary_AABB_tree.cpp
    broadphase/test_binary_AABB_tree_collision.cpp
    # geometry/shape
    geometry/shape/test_capsule.cpp
    geometry/shape/test_convex.cpp
    # geometry/shape/test_sweep_volume_geometry.cpp
    # geometry/shape/test_sweep_volume_sphere_capsule.cpp
    # geometry/shape/test_swept_volume_vs_point_collision.cpp
    # geometry/heightmap
    geometry/heightmap/test_flat_heightmap.cpp
    geometry/heightmap/test_layered_heightmap.cpp
    geometry/heightmap/test_heightmap_shape.cpp
    geometry/heightmap/test_heightmap_pair_collision.cpp
    geometry/heightmap/test_heightmap_bvh_collision.cpp
    geometry/heightmap/test_heightmap_octree2_collision.cpp
    # geometry/octree2
    geometry/octree2/test_octree_node.cpp
    geometry/octree2/test_octree_construct_by_hand.cpp
    geometry/octree2/test_octree_shape_collision.cpp
    geometry/octree2/test_octree_pair_collision.cpp
    geometry/octree2/test_octree_bvh_collision.cpp
    # narrowphase
    narrowphase/detail/test_collision_func_matrix.cpp
    narrowphase/detail/test_other_collision_result.cpp
    # narrowphase/detail/test_colliding_volume_bvh.cpp
    # narrowphase/detail/test_colliding_volume_primitive.cpp
    narrowphase/detail/primitive_shape_algorithm/test_half_space_convex.cpp
    narrowphase/detail/primitive_shape_algorithm/test_sphere_box.cpp
    narrowphase/detail/primitive_shape_algorithm/test_sphere_cylinder.cpp
    narrowphase/detail/primitive_shape_algorithm/test_fcl_box_triangle.cpp
    narrowphase/detail/primitive_shape_algorithm/test_fcl_tetrahedron.cpp
    # ccd
    ccd/test_box_pair_ccd.cpp
    ccd/test_shape_pair_interface.cpp
    ccd/test_bvh_ccd.cpp
    ccd/test_octree_shape_ccd.cpp
    ccd/test_octree_pair_ccd.cpp
    ccd/test_octree_bvh_ccd.cpp
    ccd/test_heightmap_shape_ccd.cpp
    ccd/test_heightmap_bvh_ccd.cpp
    ccd/test_heightmap_octree_ccd.cpp
    ccd/test_heightmap_pair_ccd.cpp
)

# The test that depends on internal data
if (PRIVATE_TEST_MODEL_DIR)
    list(APPEND tests private/test_fcl_point_cloud_collision.cpp)
    list(APPEND tests private/test_octree_mesh_0.cpp)
    list(APPEND tests private/test_binvox_octree.cpp)
    list(APPEND tests private/test_mesh_obj.cpp)
endif()

# Add google tests
enable_testing()
include(GoogleTest)

# Add test functor
macro(add_fcl_test test_file_name)
  # Get the name (i.e. bla.cpp => bla)
  get_filename_component(test_name ${ARGV} NAME_WE)
  add_executable(${test_name} ${ARGV})
  target_link_libraries(${test_name} PRIVATE test_fcl_utility gtest fcl)
  # add_test(${test_name} ${EXECUTABLE_OUTPUT_PATH}/${test_name})
  gtest_discover_tests(${test_name})
endmacro(add_fcl_test)

# configure location of resources
file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/fcl_resources" TEST_RESOURCES_SRC_DIR)
file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/fcl_resources" TEST_RESOURCES_BIN_DIR)
if(WIN32)
    # Correct directory separator for Windows
    string(REPLACE "\\" "\\\\" TEST_RESOURCES_SRC_DIR ${TEST_RESOURCES_SRC_DIR})
    string(REPLACE "\\" "\\\\" TEST_RESOURCES_BIN_DIR ${TEST_RESOURCES_BIN_DIR})
endif(WIN32)
configure_file("${TEST_RESOURCES_SRC_DIR}/config.h.in" "${TEST_RESOURCES_BIN_DIR}/config.h")

include_directories(.)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# The benchmark
macro(add_fcl_benchmark test_file_name)
    # Get the name (i.e. bla.cpp => bla)
    get_filename_component(test_name ${ARGV} NAME_WE)
    add_executable(${test_name} ${ARGV})
    target_link_libraries(${test_name} test_fcl_utility)
endmacro(add_fcl_benchmark)

# The general benchmark
add_fcl_benchmark(cvx_collide/gjk_benchmark.cpp)
add_fcl_benchmark(cvx_collide/mpr_benchmark.cpp)
add_fcl_benchmark(cvx_collide/mpr_refine_benchmark.cpp)
add_fcl_benchmark(geometry/heightmap/flat_heightmap_benchmark.cpp)
add_fcl_benchmark(geometry/heightmap/heightmap_shape_collision_benchmark.cpp)
add_fcl_benchmark(narrowphase/detail/primitive_shape_algorithm/benchmark_fcl_box_triangle.cpp)
add_fcl_benchmark(narrowphase/detail/primitive_shape_algorithm/benchmark_fcl_tetrahedron.cpp)

# The benchmark that depends on internal data
if (PRIVATE_TEST_MODEL_DIR)
    add_fcl_benchmark(private/octree_mesh_benchmark_0.cpp)
    add_fcl_benchmark(private/octree_pair_benchmark.cpp)
    add_fcl_benchmark(private/octree_shape_benchmark.cpp)
    add_fcl_benchmark(private/heightmap_shape_benchmark.cpp)
    add_fcl_benchmark(private/heightmap_mesh_benchmark_0.cpp)
    add_fcl_benchmark(private/sandbox_benchmark.cpp)
    add_fcl_benchmark(private/sandbox_mesh_obj.cpp)
    add_fcl_benchmark(private/octree_building_benchmark.cpp)
endif()

# Build all the tests
foreach(test ${tests})
  add_fcl_test(${test})
endforeach(test)
