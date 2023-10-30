cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(BUILD_DOC OFF)
if(NOT SUBPACKAGE)
  option(BUILD_TEST "Build tests" OFF)
  option(BUILD_PYTHON "Build python interfaces" ON)
  option(INSTALL_CMAKE_ENV "Install cmake build environment" ON)
  option(BUILD_DOC "Build documentation" OFF)
  option(BUILD_PYTHON "Build python" ON)
else(NOT SUBPACKAGE)

endif(NOT SUBPACKAGE)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_DOUBLE "Build double support" ON)
option(BUILD_INT "Build int support" ON)
option(BUILD_COMPLEX "Build complex support" ON)
option(BUILD_BYTE "Build byte support" ON)
option(BUILD_SHORT "Build short support" ON)

if(BUILD_TEST)
	enable_testing()
endif(BUILD_TEST)

if(BUILD_DOUBLE)
  set(USE_DOUBLE true)
endif(BUILD_DOUBLE)

if(BUILD_COMPLEX)
  set(USE_COMPLEX true)
endif(BUILD_COMPLEX)

if(BUILD_INT)
  set(USE_INT true)
endif(BUILD_INT)

if(BUILD_BYTE)
  set(USE_BYTE true)
endif(BUILD_BYTE)
if(BUILD_SHORT)
  set(USE_SHORT true)
endif(BUILD_SHORT)

if(POLICY CMP0028)
  cmake_policy(SET CMP0028 NEW)
endif()

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

project(sepVector LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE True)

project(sepVector)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

# configure a header file to pass some of the CMake settings
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Boost REQUIRED)
find_package(TBB REQUIRED)
set(sepVector_MAJOR_VERSION 0)
set(sepVector_MINOR_VERSION 2)
set(sepVector_VERSION_STRING
    ${sepVector_MAJOR_VERSION}.${sepVector_MINOR_VERSION})

set(sepVector_DIRNAME_lib "lib")
set(sepVector_DIRNAME_include "include")

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

if(NOT SUBPACKAGE)
  execute_process(COMMAND git submodule update --init --recursive --
                          external/python-solver
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  execute_process(COMMAND git submodule update --init --recursive --
                          external/hypercube
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

configure_file("${PROJECT_SOURCE_DIR}/sepVectorConfig.h.in"
               "${PROJECT_BINARY_DIR}/sepVectorConfig.h")

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
    add_custom_target(doc_doxygen_sepVector ALL
                      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)
  else(DOXYGEN_FOUND)
    message(
      "Doxygen need to be installed to generate the doxygen documentation")
  endif(DOXYGEN_FOUND)
else(BUILD_DOC)
endif(BUILD_DOC)
set(SUBPACKAGE true)
if(NOT TARGET hypercube)
add_subdirectory(external/hypercube)
endif(NOT TARGET hypercube)
add_subdirectory(external/python-solver/GenericSolver)
add_subdirectory(lib)
if(BUILD_PYTHON)
  add_subdirectory(python)
endif(BUILD_PYTHON)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  sepVectorConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/sepVectorConfig.cmake"
  INSTALL_DESTINATION
  cmake/SEP/IO)

install(EXPORT sepVectorTargets
        FILE sepVectorTargets.cmake
        NAMESPACE SEP::IO
        DESTINATION cmake/sep/io)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${sepVector_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${sepVector_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_NAME "sepVector")
set(CPACK_PACKAGE_CONTACT "Robert Clapp")
include(CPack)
