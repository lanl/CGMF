# copy reference files to local
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/yields.cgmf.serial.0.reference
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


# cgmf.x arguments
set(args     -n -1000 -e 10.0 -i 92235 -f yields.cgmf.serial)
set(args.mpi -n -1000 -e 10.0 -i 92235 -f yields.cgmf.parallel)


if(NOT ${util}.x.MPI)
  # cgmf.x serial tests
  add_test(NAME cgmf_u235nf-10-yields_init
           COMMAND ${CMAKE_COMMAND} -E remove
                   yields.cgmf.serial.0)
  set_tests_properties(cgmf_u235nf-10-yields_init PROPERTIES
                       FIXTURES_SETUP cgmf_u235nf-10-yields_init)

  add_test(NAME cgmf_u235nf-10-yields_exec
           COMMAND ${util}.x ${args}
           COMMAND_EXPAND_LISTS)
  set_tests_properties(cgmf_u235nf-10-yields_exec PROPERTIES
                       FIXTURES_SETUP cgmf_u235nf-10-yields_exec
                       FIXTURES_REQUIRED cgmf_u235nf-10-yields_init)

  add_test(NAME cgmf_u235nf-10-yields_diff
           COMMAND ${CMAKE_COMMAND} -E compare_files
                   yields.cgmf.serial.0.reference yields.cgmf.serial.0)
  set_tests_properties(cgmf_u235nf-10-yields_diff PROPERTIES
                       FIXTURES_REQUIRED cgmf_u235nf-10-yields_exec)

else(${util}.x.MPI)
  # cgmf.mpi.x mpi tests
  set(numprocs.mpi 2)
  add_test(NAME cgmf_u235nf-10-yields_mpi-init
           COMMAND ${CMAKE_COMMAND} -E remove
                   yields.cgmf.parallel.0)
  set_tests_properties(cgmf_u235nf-10-yields_mpi-init PROPERTIES
                       FIXTURES_SETUP cgmf_u235nf-10-yields_mpi-init)


  add_test(NAME cgmf_u235nf-10-yields_mpi-exec
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                   ${numprocs.mpi} ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:${util}.x> ${MPIEXEC_POSTFLAGS} ${args.mpi}
           COMMAND_EXPAND_LISTS)
  set_tests_properties(cgmf_u235nf-10-yields_mpi-exec PROPERTIES
                       FIXTURES_SETUP cgmf_u235nf-10-yields_mpi-exec
                       FIXTURES_REQUIRED cgmf_u235nf-10-yields_mpi-init)

  add_test(NAME cgmf_u235nf-10-yields_mpi-diff-0
           COMMAND ${CMAKE_COMMAND} -E compare_files
                   yields.cgmf.serial.0.reference yields.cgmf.parallel.0)
  set_tests_properties(cgmf_u235nf-10-yields_mpi-diff-0 PROPERTIES
                       FIXTURES_REQUIRED cgmf_u235nf-10-yields_mpi-exec)

endif()
