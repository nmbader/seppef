cmake_minimum_required(VERSION 3.5 FATAL_ERROR)


if(NOT SUBPACKAGE)
  option(BUILD_PYTHON "Build python interfaces" ON)
  option(INSTALL_CMAKE_ENV "Install cmake build environment" ON)
  option(BUILD_DOC "Build documentation" OFF)
  option(BUILD_TEST "Build with test support" OFF)
else(NOT SUBPACKAGE)

endif(NOT SUBPACKAGE)

project(hypercube LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE True)
# configure a header file to pass some of the CMake settings

set(hypercube_MAJOR_VERSION 0)
set(hypercube_MINOR_VERSION 2)
set(hypercube_VERSION_STRING
    ${hypercube_MAJOR_VERSION}.${hypercube_MINOR_VERSION})

include_directories("${PROJECT_BINARY_DIR}")

set(hypercube_DIRNAME_cmake "lib/cmake")
set(hypercube_DIRNAME_lib "lib")
set(hypercube_DIRNAME_include "include")

include(CMakePackageConfigHelpers)

# first we can indicate the documentation build as an option and set it to ON by
# default

if(BUILD_DOC)
  # check if Doxygen is installed
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # note the option ALL which allows to build the docs together with the
    # application
    add_custom_target(doc_doxygen ALL
                      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)
  else(DOXYGEN_FOUND)
    message(
      "Doxygen need to be installed to generate the doxygen documentation")
  endif(DOXYGEN_FOUND)
endif(BUILD_DOC)

if(NOT SUBPACKAGE)
  execute_process(COMMAND git submodule update --init --recursive --
                          external/pybind11
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  if(NOT BUILD_GCP)
    execute_process(COMMAND git submodule update --init --recursive --
                            external/googletest
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif(NOT BUILD_GCP)
endif(NOT SUBPACKAGE)

set(SUBPACKAGE true)

if(BUILD_TEST)
  enable_testing()
endif(BUILD_TEST)

if(NOT BUILD_GCP AND BUILD_TEST AND NOT TARGET gtest_main)

  add_subdirectory(external/googletest)
endif()
if(NOT TARGET pybind11)
add_subdirectory(external/pybind11)
else()
endif(NOT TARGET pybind11)
add_subdirectory(lib)
add_subdirectory(python)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  hypercubeConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/hypercubeConfig.cmake"
  INSTALL_DESTINATION
  cmake/SEP/IO)

install(EXPORT HypercubeTargets
        FILE hypercubeTargets.cmake
        NAMESPACE SEP::IO
        DESTINATION cmake/sep/io)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${hypercube_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${hypercube_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_NAME "hypercube")
set(CPACK_PACKAGE_CONTACT "Robert Clapp")
include(CPack)
