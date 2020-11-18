set(util_src ${CMAKE_CURRENT_SOURCE_DIR}/cgmf.cpp)

get_filename_component(util ${util_src} NAME_WE)

add_executable(${util}.x ${util_src})
set_target_properties(${util}.x PROPERTIES COMPILE_FLAGS "${CGMF_CXX_FLAGS}")
set_target_properties(${util}.x PROPERTIES LINK_FLAGS "${CGMF_EXE_LINK_FLAGS}")
set_target_properties(${util}.x PROPERTIES COMPILE_DEFINITIONS "${CGMF_COMPILE_DEFS}")
target_link_libraries(${util}.x cgmf)

# Serial or MPI Build
option(${util}.x.MPI "Build CGMF executable with MPI parallelism" ${MPI})
if(${util}.x.MPI)
  set(MPI_HOME "" CACHE STRING "Specify MPI_HOME for specific MPI implementation ( MPI_ROOT )" )
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  target_compile_definitions(${util}.x PUBLIC MPIRUN)
  target_link_libraries(${util}.x ${MPI_CXX_LIBRARIES})
  set_target_properties(${util}.x PROPERTIES OUTPUT_NAME ${util}.mpi.x)
endif()

if(cgmf.tests)
  add_subdirectory(tests)
endif()

install(TARGETS ${util}.x DESTINATION bin)
