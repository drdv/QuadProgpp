cmake_minimum_required(VERSION 2.8.12)
project(QuadProg++ CXX)

set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

option(QUADPROGPP_BUILD_EXAMPLE     "Build example"     ON)
option(QUADPROGPP_ENABLE_TRACING    "Enable tracing"    OFF)
option(QUADPROGPP_ENABLE_EIGEN      "Enable Eigen"      OFF)

#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_BUILD_TYPE Release)
set (CMAKE_VERBOSE_MAKEFILE     ON)

if(QUADPROGPP_ENABLE_EIGEN)
    # Search for Eigen
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif(QUADPROGPP_ENABLE_EIGEN)


set(QUADPROGPP_SOURCE_DIR   "${PROJECT_SOURCE_DIR}/src")

configure_file("cmake/config.hh.in"        "${QUADPROGPP_SOURCE_DIR}/config.hh")

add_subdirectory(${QUADPROGPP_SOURCE_DIR})
