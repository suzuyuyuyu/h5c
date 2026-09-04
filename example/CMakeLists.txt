cmake_minimum_required(VERSION 3.20)
project(h5c_examples LANGUAGES C)

# The examples are a STANDALONE project on purpose: the library build never
# references them. Consuming h5c through find_package() here is the same thing
# a real user writes, so the examples cannot drift away from reality.
#
#   cmake -S example -B example/build/intel \
#         -DCMAKE_PREFIX_PATH="$HOME/.local/opt/intel/h5c-0.1.0"
#   cmake --build example/build/intel
#
# When h5c is not installed, fall back to the sibling source tree so the
# examples still build in a fresh checkout.
find_package(h5c CONFIG QUIET)
if(NOT h5c_FOUND)
    message(STATUS "h5c not installed; building it from ${CMAKE_CURRENT_SOURCE_DIR}/..")
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(H5C_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    add_subdirectory(".." "${CMAKE_CURRENT_BINARY_DIR}/h5c")
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(h5c_example_targets)

function(h5c_add_example target source)
    add_executable("${target}" "${source}")
    target_link_libraries("${target}" PRIVATE h5c::h5c)
    list(APPEND h5c_example_targets "${target}")
    set(h5c_example_targets "${h5c_example_targets}" PARENT_SCOPE)
endfunction()

h5c_add_example(example_serial      serial/main.c)
h5c_add_example(example_serial_interleaved interleaved/main.c)

# The parallel example needs an h5c built with H5C_ENABLE_PARALLEL=ON, which
# puts H5C_HAVE_PARALLEL on the target. Checking the target rather than a path
# works for both an installed h5c and an add_subdirectory() one.
get_target_property(_h5c_defs h5c::h5c INTERFACE_COMPILE_DEFINITIONS)
if(_h5c_defs AND "H5C_HAVE_PARALLEL" IN_LIST _h5c_defs)
    find_package(MPI REQUIRED COMPONENTS C)
    h5c_add_example(example_parallel parallel/main.c)
    target_link_libraries(example_parallel PRIVATE MPI::MPI_C)
    h5c_add_example(example_parallel_interleaved parallel-interleaved/main.c)
    target_link_libraries(example_parallel_interleaved PRIVATE MPI::MPI_C)
    h5c_add_example(example_parallel_visualization visualization/main.c)
    target_link_libraries(example_parallel_visualization PRIVATE MPI::MPI_C m)
else()
    message(STATUS "h5c has no parallel support; skipping example_parallel")
endif()

add_custom_target(examples DEPENDS ${h5c_example_targets})
