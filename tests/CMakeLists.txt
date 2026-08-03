if (EC_INIT_DONE LESS 2)
  include(CTest)
  add_custom_target(check COMMAND CTEST_OUTPUT_ON_FAILURE=TRUE ${CMAKE_CTEST_COMMAND} -V)
  add_custom_target(test_execs)
  add_dependencies(check test_execs)

  if(BUILD_TESTING)
    message(STATUS "(EC) Generating random_tools tests")

    add_executable(c_demo EXCLUDE_FROM_ALL random_demo.c)
    target_link_libraries(c_demo random m)
    target_include_directories(c_demo PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
    add_dependencies(c_demo random)
    add_dependencies(test_execs c_demo)
    add_test(NAME c_demo COMMAND ${CMAKE_CURRENT_BINARY_DIR}/c_demo)

    add_executable(f_demo EXCLUDE_FROM_ALL random_demo.F90)
    find_package(OpenMP)
    target_link_libraries(f_demo
      $<$<LINK_LANGUAGE:Fortran>:OpenMP::OpenMP_Fortran>
      $<$<LINK_LANGUAGE:C>:OpenMP::OpenMP_C>)
    target_link_libraries(f_demo random m)
    target_include_directories(f_demo PRIVATE ${CMAKE_BINARY_DIR}/src)
    add_dependencies(f_demo random randomfunctionsINC)
    add_dependencies(test_execs f_demo)
    add_test(NAME f_demo COMMAND ${CMAKE_CURRENT_BINARY_DIR}/f_demo)

    if (MPI_FOUND)
      # Check MPI implementation to be able to give the proper flags when running tests
      find_program(RPN_MPIRUN mpirun)
      if (${RPN_MPIRUN} STREQUAL "RPN_MPIRUN-NOTFOUND")
        message(FATAL_ERROR "mpirun not found. Won't be able to configure the version.")
      else()
        execute_process(
          COMMAND ${RPN_MPIRUN} --version 2>&1
          COMMAND grep "Open MPI"
          COMMAND grep -o "[0-9]\\.[0-9]"
          OUTPUT_VARIABLE MPI_IMPLEMENTATION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (MPI_IMPLEMENTATION MATCHES "[0-9]\\.[0-9]")
          set(EXTRA_MPI_FLAGS ${EXTRA_MPI_FLAGS} --oversubscribe)
          if (CMAKE_Fortran_COMPILER_ID MATCHES ".*Flang")
            set(EXTRA_MPI_FLAGS ${EXTRA_MPI_FLAGS} -x LD_LIBRARY_PATH)
          else()
            set(EXTRA_MPI_FLAGS ${EXTRA_MPI_FLAGS} --bind-to none)
          endif()
        endif()
      endif()
      
      foreach(GENERIC_TEST R250 MT19937 SHR3 XSR128 XSR128R)
	add_executable(random_${GENERIC_TEST} EXCLUDE_FROM_ALL random_generic_test.c)
	target_compile_definitions(random_${GENERIC_TEST} PRIVATE TEST_${GENERIC_TEST})
	target_link_libraries(random_${GENERIC_TEST} PRIVATE random m MPI::MPI_C)
	target_include_directories(random_${GENERIC_TEST} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
	add_dependencies(random_${GENERIC_TEST} random)
	add_dependencies(test_execs random_${GENERIC_TEST})
	add_test(NAME random_${GENERIC_TEST} COMMAND ${MPIEXEC_EXECUTABLE} ${EXTRA_MPI_FLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/random_${GENERIC_TEST})
      endforeach()

      add_executable(random_gaussian EXCLUDE_FROM_ALL random_gaussian_test.c)
      target_compile_definitions(random_gaussian PUBLIC FULL_TEST)
      target_link_libraries(random_gaussian random m MPI::MPI_C)
      target_include_directories(random_gaussian PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
      add_dependencies(random_gaussian random)
      add_dependencies(test_execs random_gaussian)
      add_test(NAME random_gaussian COMMAND ${MPIEXEC_EXECUTABLE} ${EXTRA_MPI_FLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/random_gaussian)

      add_executable(random_gaussian_profile EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/src/random_gaussian.c random_gaussian_test.c)
      target_compile_definitions(random_gaussian_profile PUBLIC FULL_TEST PROFILE)
      target_link_libraries(random_gaussian_profile random m MPI::MPI_C)
      target_include_directories(random_gaussian_profile PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
      add_dependencies(random_gaussian_profile random)
      add_dependencies(test_execs random_gaussian_profile)
      add_test(NAME random_gaussian_profile COMMAND ${MPIEXEC_EXECUTABLE} ${EXTRA_MPI_FLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/random_gaussian_profile)
    endif()
  endif()
endif()
