# -*- mode: cmake -*-
# vi: set ft=cmake :

# Copyright (c) 2012, Willow Garage, Inc.
# Copyright (c) 2016, Toyota Research Institute, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.10)

project(fcl CXX)

# Compiler ID for Apple Clang is now AppleClang.
if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif()

option(BUILD_HEADER_ONLY              "As head-only lib"                               OFF)
option(FCL_STATIC_LIBRARY             "If not built as header-only, static/shared"     ON )
option(BUILD_TESTING                  "Build FCL Testing"                              ON )
option(FCL_TREAT_WARNINGS_AS_ERRORS   "Treat warnings as errors"                       OFF)

# set the default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build; options are Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# This shouldn't be necessary, but there has been trouble
# with MSVC being set off, but MSVCXX ON.
if(MSVC OR MSVC90 OR MSVC10)
    set(MSVC ON)
endif (MSVC OR MSVC90 OR MSVC10)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)
include(GNUInstallDirs)
include(CompilerSettings)
include(FCLVersion)

#===============================================================================
# Detect SSE flags
#===============================================================================
set(SSE_FLAGS "")
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(SSE_FLAGS -mfpmath=sse -msse -msse2 -msse3 -mssse3)
elseif(MSVC)
  # Win64 will add the flag automatically
  if(CMAKE_VS_PLATFORM_NAME STREQUAL "Win32")
    set(SSE_FLAGS /arch:SSE2)
  endif()
endif()

# Whether to enable SSE
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set(FCL_TARGET_SUPPORT_X64_SSE OFF)
else()
  # Does the host support SSE
  cmake_host_system_information(RESULT _has_sse QUERY HAS_SSE)
  cmake_host_system_information(RESULT _has_sse2 QUERY HAS_SSE2)

  # Does the compiler support SSE
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("${SSE_FLAGS}" _compiler_supports_sse)

  if(_has_sse AND _has_sse2 AND _compiler_supports_sse)
    set(FCL_TARGET_SUPPORT_X64_SSE ON)
  else()
    set(FCL_TARGET_SUPPORT_X64_SSE OFF)
  endif()
endif()


option(FCL_USE_X64_SSE "Whether FCL should x64 SSE instructions" ${FCL_TARGET_SUPPORT_X64_SSE})
if(FCL_USE_X64_SSE)
  message(STATUS "FCL uses SSE")
  add_compile_options(${SSE_FLAGS})
else()
  # Unset SSE_FLAGS so it doesn't get used in subdirectories
  message(STATUS "FCL does not use SSE")
  set(SSE_FLAGS "")
endif()

option(FCL_USE_HOST_NATIVE_ARCH "Whether FCL should use cflags from the host used to compile" OFF)
if (FCL_USE_HOST_NATIVE_ARCH)
  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  else()
    message(WARNING "FCL_USE_HOST_NATIVE_ARCH is only supported in Linux. No effect.")
  endif()
endif()

#===============================================================================
# Find required dependency Eigen3 (>= 3.0.5)
#
# If Eigen3 is not found, manually set the cache variable EIGEN3_INCLUDE_DIR
#===============================================================================
find_package(Eigen3 3.0.5 QUIET CONFIG)

# If Eigen3Config.cmake is not found, use the FindEigen3.cmake module
if(NOT Eigen3_FOUND)
  find_package(Eigen3 3.0.5 QUIET MODULE)
  set(Eigen3_FOUND ON)
endif()

if(Eigen3_FOUND)
  set(FCL_HAVE_EIGEN TRUE)
else()
  message(SEND_ERROR "EIGEN3 (>= 3.0.5) is required by FCL")
  set(FCL_HAVE_EIGEN FALSE)
endif()

# --------------------------------------------------------------------
# Install/uninstall targets
# --------------------------------------------------------------------
if(TARGET Eigen3::Eigen)
  set(FIND_DEPENDENCY_EIGEN3 "find_dependency(Eigen3)")
else()
  set(FIND_DEPENDENCY_EIGEN3)
endif()

if(WIN32 AND NOT CYGWIN)
  set(FCL_INSTALL_CONFIGDIR CMake)
else()
  set(FCL_INSTALL_CONFIGDIR share/cmake/${PROJECT_NAME})
endif()

configure_package_config_file(fcl-config.cmake.in ${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${FCL_INSTALL_CONFIGDIR}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR
)

write_basic_package_version_file(${PROJECT_NAME}-config-version.cmake
  VERSION ${FCL_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION ${FCL_INSTALL_CONFIGDIR}
  COMPONENT Development
)

set(INSTALL_EXTERN_HEADER_FOR_LIB FALSE)
if (NOT BUILD_HEADER_ONLY)
    set(INSTALL_EXTERN_HEADER_FOR_LIB TRUE)
endif()

# FCL's own include dir should be at the front of the include path
include_directories(BEFORE "include")
add_subdirectory(include/fcl)

# Define fcl lib (header-only)
if (BUILD_HEADER_ONLY)
  message(STATUS "Use fcl as a header-only lib. Might increase the build time in your application")
  add_library(${PROJECT_NAME} INTERFACE)
  target_include_directories(${PROJECT_NAME}
          INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  # Be sure to pass to the consumer the set of SIMD used in the compilation
  target_compile_options(${PROJECT_NAME} INTERFACE ${SSE_FLAGS})

  # Use the IMPORTED target from newer versions of Eigen3Config.cmake if
  # available, otherwise fall back to EIGEN3_INCLUDE_DIRS from older versions of
  # Eigen3Config.cmake or EIGEN3_INCLUDE_DIR from FindEigen3.cmake
  if(TARGET Eigen3::Eigen)
    # Note that Eigen3::Eigen is an INTERFACE library, so the INCLUDE_DIRECTORIES
    # and INTERFACE_INCLUDE_DIRECTORIES are populated, but nothing is actually
    # linked
    target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
  elseif(EIGEN3_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME} INTERFACE "${EIGEN3_INCLUDE_DIRS}")
  else()
    target_include_directories(${PROJECT_NAME} INTERFACE "${EIGEN3_INCLUDE_DIR}")
  endif()
else ()
  message(STATUS "Build the fcl lib.")
  add_subdirectory(src)
endif ()

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT Runtime ${NAMELINK_COMPONENT_OPTION}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT Runtime
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

export(EXPORT ${PROJECT_NAME}-targets
        NAMESPACE fcl::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake")

export(PACKAGE ${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}-targets
        NAMESPACE fcl::
        DESTINATION "${FCL_INSTALL_CONFIGDIR}")

# The tests
if(BUILD_TESTING)
  add_subdirectory(test)
endif()