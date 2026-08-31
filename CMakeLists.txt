cmake_minimum_required(VERSION 3.20)

project(h5c
    VERSION 0.1.0
    DESCRIPTION "Compact C wrapper around the HDF5 C API"
    LANGUAGES C
)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(H5C_ENABLE_PARALLEL "Build the parallel (MPI) API" OFF)
option(H5C_BUILD_TESTS "Build the test programs" ${BUILD_TESTING})

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

# ---------------------------------------------------------------------
# library
# ---------------------------------------------------------------------

set(h5c_sources
    src/h5c_error.c
    src/h5c_type.c
    src/h5c_file.c
    src/h5c_dataset.c
    src/h5c_string.c
    src/h5c_attribute.c
    src/h5c_interleaved.c
)
if(H5C_ENABLE_PARALLEL)
    list(APPEND h5c_sources src/h5c_parallel.c)
endif()

add_library(h5c ${h5c_sources})
add_library(h5c::h5c ALIAS h5c)

target_include_directories(h5c
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# HDF5 requirements are carried on the target so dependents inherit them.
if(TARGET HDF5::HDF5)
    target_link_libraries(h5c PUBLIC HDF5::HDF5)
else()
    target_include_directories(h5c PUBLIC ${HDF5_C_INCLUDE_DIRS})
    target_link_libraries(h5c PUBLIC ${HDF5_C_LIBRARIES})
endif()

if(H5C_ENABLE_PARALLEL)
    target_link_libraries(h5c PUBLIC MPI::MPI_C)
    target_compile_definitions(h5c PUBLIC H5C_HAVE_PARALLEL)
endif()

set_target_properties(h5c PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# ---------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------

if(H5C_BUILD_TESTS)
    add_subdirectory(test)
endif()

# ---------------------------------------------------------------------
# install / export
# ---------------------------------------------------------------------

install(TARGETS h5c
    EXPORT h5cTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
