cmake_minimum_required(VERSION 2.8.12)
project(QuadProg++ CXX)

option(QUADPROGPP_BUILD_EXAMPLE     "Build example" ON)
option(QUADPROGPP_ENABLE_TRACING    "Enable tracing" OFF)


set(QUADPROGPP_SOURCE_DIR   "${PROJECT_SOURCE_DIR}/src")

configure_file("config.hh.in"        "${QUADPROGPP_SOURCE_DIR}/config.hh")

add_subdirectory(${QUADPROGPP_SOURCE_DIR})
