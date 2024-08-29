file(GLOB_RECURSE FCL_SOURCE_CODE ${CONFIGURE_DEPENDS_OPTION}
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

if(FCL_STATIC_LIBRARY)
    add_library(${PROJECT_NAME} STATIC ${FCL_HEADERS} ${FCL_SOURCE_CODE})
else()
    add_library(${PROJECT_NAME} SHARED ${FCL_HEADERS} ${FCL_SOURCE_CODE})
endif()

target_include_directories(${PROJECT_NAME}
        INTERFACE
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Be sure to pass to the consumer the set of SIMD used in the compilation
target_compile_options(${PROJECT_NAME} PUBLIC ${SSE_FLAGS})

# Use the IMPORTED target from newer versions of Eigen3Config.cmake if
# available, otherwise fall back to EIGEN3_INCLUDE_DIRS from older versions of
# Eigen3Config.cmake or EIGEN3_INCLUDE_DIR from FindEigen3.cmake
if(TARGET Eigen3::Eigen)
    # Note that Eigen3::Eigen is an INTERFACE library, so the INCLUDE_DIRECTORIES
    # and INTERFACE_INCLUDE_DIRECTORIES are populated, but nothing is actually
    # linked
    target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
elseif(EIGEN3_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME} PUBLIC "${EIGEN3_INCLUDE_DIRS}")
else()
    target_include_directories(${PROJECT_NAME} PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()