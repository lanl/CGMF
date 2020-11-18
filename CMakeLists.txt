cmake_minimum_required(VERSION 3.16.2)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Specifies the build type on single-configuration generators.")
endif()

if(NOT DEFINED PROJECT_NAME)
  set(NOT_SUBPROJECT ON)
endif()

file(STRINGS VERSION CGMF_VER)
set(cgmf_version ${CGMF_VER} CACHE STRING "CGMF Library Version")

project(cgmf VERSION ${CGMF_VER} LANGUAGES C CXX)

# Set C++11 Standard
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(libcgmf)

option(cgmf.tests "Configure tests for CGMF" ON)

if(cgmf.tests)
  enable_testing()
endif()

if(NOT_SUBPROJECT)
  add_subdirectory(utils)
endif()

install(DIRECTORY data DESTINATION cgmf)
