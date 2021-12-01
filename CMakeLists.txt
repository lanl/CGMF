cmake_minimum_required(VERSION 3.16.2)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Specifies the build type on single-configuration generators.")
endif()

if(DEFINED PROJECT_NAME)
  set(SUBPROJECT ON)
  set(INSTALL_SUBPROJECTS ON CACHE BOOL "Install subproject dependencies")
else()
  set(SUBPROJECT OFF)
  set(INSTALL_SUBPROJECTS ON CACHE BOOL "Install subproject dependencies" FORCE)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

file(STRINGS VERSION CGMF_VER)
set(cgmf_version ${CGMF_VER} CACHE STRING "CGMF Library Version")

project(cgmf VERSION ${CGMF_VER} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)

include(CMakeDependentOption)
include(GNUInstallDirs)

CMAKE_DEPENDENT_OPTION(cgmf.installation 
                       "Install CGMF (default is ON)" ON
                       "NOT SUBPROJECT" ${INSTALL_SUBPROJECTS}
)
mark_as_advanced(cgmf.installation)

add_subdirectory(libcgmf)

option(cgmf.tests "Configure tests for CGMF" ON)

if(cgmf.tests)
  enable_testing()
endif()

if(NOT SUBPROJECT)
  add_subdirectory(utils)
endif()

if(cgmf.installation)
  install(DIRECTORY data/
          DESTINATION "${CMAKE_INSTALL_DATADIR}/data/cgmf-${CGMF_VER}"
          FILES_MATCHING PATTERN "*.dat*"
  )
endif()
