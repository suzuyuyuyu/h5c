cmake_minimum_required(VERSION 3.20)

project(h5c
    VERSION 2.0.0
    DESCRIPTION "Compact C wrapper around the HDF5 C API"
    LANGUAGES C
)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(H5C_ENABLE_PARALLEL "Build the parallel (MPI) API" OFF)
option(H5C_BUILD_TESTS "Build the test programs" ${BUILD_TESTING})

# Layout version of the visualization HDF5 this build writes. Shared with
# h5fortran as a format contract, so it is set here rather than derived from
# PROJECT_VERSION -- see include/h5c/h5c_version.h.in for why.
set(H5C_SCHEME_VERSION 1)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/h5c/h5c_version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/h5c/h5c_version.h"
    @ONLY)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(H5C_ENABLE_PARALLEL)
    set(HDF5_PREFER_PARALLEL ON)
endif()
find_package(HDF5 REQUIRED COMPONENTS C)

if(H5C_ENABLE_PARALLEL)
    find_package(MPI REQUIRED COMPONENTS C)
    if(DEFINED HDF5_IS_PARALLEL AND NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR
            "H5C_ENABLE_PARALLEL requires an HDF5 build with parallel I/O support")
    endif()
endif()

set(h5c_serial_sources
    src/serial/h5c_error.c
    src/serial/h5c_type.c
    src/serial/h5c_file.c
    src/serial/h5c_dataset.c
    src/serial/h5c_string.c
    src/serial/h5c_attribute.c
    src/serial/h5c_interleaved.c
    src/serial/h5c_viz.c
)
add_library(h5c_serial ${h5c_serial_sources})
add_library(h5c::h5c_serial ALIAS h5c_serial)
add_library(h5c INTERFACE)
add_library(h5c::h5c ALIAS h5c)
target_link_libraries(h5c INTERFACE h5c_serial)
set(h5c_install_targets h5c h5c_serial)

target_include_directories(h5c_serial
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# HDF5 requirements are carried on the target so dependents inherit them.
if(TARGET HDF5::HDF5)
    target_link_libraries(h5c_serial PUBLIC HDF5::HDF5)
else()
    target_include_directories(h5c_serial PUBLIC ${HDF5_C_INCLUDE_DIRS})
    target_link_libraries(h5c_serial PUBLIC ${HDF5_C_LIBRARIES})
endif()

if(H5C_ENABLE_PARALLEL)
    add_library(h5c_parallel src/parallel/parallel.c src/parallel/viz_mpi.c)
    add_library(h5c::h5c_parallel ALIAS h5c_parallel)
    target_include_directories(h5c_parallel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(h5c_parallel PUBLIC h5c_serial MPI::MPI_C)
    target_compile_definitions(h5c_parallel PUBLIC H5C_HAVE_PARALLEL)
    target_link_libraries(h5c INTERFACE h5c_parallel)
    list(APPEND h5c_install_targets h5c_parallel)
endif()

set_target_properties(h5c_serial PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)
if(H5C_ENABLE_PARALLEL)
    set_target_properties(h5c_parallel PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
endif()

if(H5C_BUILD_TESTS)
    add_subdirectory(test)
endif()

install(TARGETS ${h5c_install_targets}
    EXPORT h5cTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PATTERN "*.in" EXCLUDE)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/h5c/h5c_version.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/h5c)

install(EXPORT h5cTargets
    FILE h5cTargets.cmake
    NAMESPACE h5c::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/h5c
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/h5cConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/h5cConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/h5c
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/h5cConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/h5cConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/h5cConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/h5c
)
