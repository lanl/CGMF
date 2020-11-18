# copy reference files to local
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/histories.cgmf.parallel.0.reference
          ${CMAKE_CURRENT_SOURCE_DIR}/histories.cgmf.parallel.1.reference
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(READ ${CMAKE_CURRENT_BINARY_DIR}/histories.cgmf.parallel.0.reference CONTENTS)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/histories.cgmf.serial.0.reference "${CONTENTS}")
file(READ ${CMAKE_CURRENT_BINARY_DIR}/histories.cgmf.parallel.1.reference CONTENTS)
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/histories.cgmf.serial.0.reference "${CONTENTS}")


# cgmf.x arguments
set(args     -n 40 -e 10.0 -i 92238 -f histories.cgmf.serial)
set(args.mpi -n 20 -e 10.0 -i 92238 -f histories.cgmf.parallel)


if(NOT ${util}.x.MPI)
  # cgmf.x serial tests
  add_test(NAME cgmf_u238nf-10-events_init
           COMMAND ${CMAKE_COMMAND} -E remove
                   histories.cgmf.serial.0)
  set_tests_properties(cgmf_u238nf-10-events_init PROPERTIES
                       FIXTURES_SETUP cgmf_u238nf-10-events_init)

  add_test(NAME cgmf_u238nf-10-events_exec
           COMMAND ${util}.x ${args}
           COMMAND_EXPAND_LISTS)
  set_tests_properties(cgmf_u238nf-10-events_exec PROPERTIES
                       FIXTURES_SETUP cgmf_u238nf-10-events_exec
                       FIXTURES_REQUIRED cgmf_u238nf-10-events_init)

  add_test(NAME cgmf_u238nf-10-events_diff
           COMMAND ${CMAKE_COMMAND} -E compare_files
                   histories.cgmf.serial.0.reference histories.cgmf.serial.0)
  set_tests_properties(cgmf_u238nf-10-events_diff PROPERTIES
                       FIXTURES_REQUIRED cgmf_u238nf-10-events_exec)

else(${util}.x.MPI)
  # cgmf.mpi.x mpi tests
  set(numprocs.mpi 2)
  add_test(NAME cgmf_u238nf-10-events_mpi-init
           COMMAND ${CMAKE_COMMAND} -E remove
                   histories.cgmf.parallel.0 histories.cgmf.parallel.1)
  set_tests_properties(cgmf_u238nf-10-events_mpi-init PROPERTIES
                       FIXTURES_SETUP cgmf_u238nf-10-events_mpi-init)


  add_test(NAME cgmf_u238nf-10-events_mpi-exec
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                   ${numprocs.mpi} ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:${util}.x> ${MPIEXEC_POSTFLAGS} ${args.mpi}
           COMMAND_EXPAND_LISTS)
  set_tests_properties(cgmf_u238nf-10-events_mpi-exec PROPERTIES
                       FIXTURES_SETUP cgmf_u238nf-10-events_mpi-exec
                       FIXTURES_REQUIRED cgmf_u238nf-10-events_mpi-init)

  add_test(NAME cgmf_u238nf-10-events_mpi-diff-0
           COMMAND ${CMAKE_COMMAND} -E compare_files
                   histories.cgmf.parallel.0.reference histories.cgmf.parallel.0)
  set_tests_properties(cgmf_u238nf-10-events_mpi-diff-0 PROPERTIES
                       FIXTURES_REQUIRED cgmf_u238nf-10-events_mpi-exec)

  add_test(NAME cgmf_u238nf-10-events_mpi-diff-1
           COMMAND ${CMAKE_COMMAND} -E compare_files
                   histories.cgmf.parallel.1.reference histories.cgmf.parallel.1)
  set_tests_properties(cgmf_u238nf-10-events_mpi-diff-1 PROPERTIES
                       FIXTURES_REQUIRED cgmf_u238nf-10-events_mpi-exec)

endif()
