cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(genericIO LANGUAGES CXX)


set(GEN_IO_TYPES \"JSON\",\"FILE_BUFFERS\")

if(POLICY CMP0028)
  cmake_policy(SET CMP0028 NEW)
endif()

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

if(NOT SUBPACKAGE)
  option(BUILD_DOC "Build documentation" OFF)
  option(BUILD_PYTHON "Build python interfaces" ON)
  option(INSTALL_CMAKE_ENV "Install cmake build environment" ON)
  option(BUILD_UTIL "Build utilities" ON)
  option(BUILD_TEST "Build with test support" OFF)
  option(BUILD_DEBUG_LIBS "Build gperftool libs" OFF)
  option(BUILD_GCP "Build GCP" OFF)
  option(BUILD_SEP "Build with SEPlib support if found" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)

else(NOT SUBPACKAGE)
  option(BUILD_UTIL "Build utilities" ON)
endif(NOT SUBPACKAGE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE True)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(TBB REQUIRED)

set(genericIO_MAJOR_VERSION 0)
set(genericIO_MINOR_VERSION 4)
set(genericIO_VERSION_STRING
    ${genericIO_MAJOR_VERSION}.${genericIO_MINOR_VERSION})

set(genericIO_VERSION ${genericIO_VERSION_STRING})
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

if(NOT SUBPACKAGE)

if(NOT TARGET sepVector)
  execute_process(COMMAND git submodule update --init --recursive --
                          external/sepVector
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif(NOT TARGET sepVector)

  if(NOT TARGET buffers)
  execute_process(COMMAND git submodule update --init --recursive --
                          external/buffers
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif(NOT TARGET buffers)
else()
endif(NOT SUBPACKAGE)

find_path(RSFINCDIR rsf.h ${RSF_DIR}/include)

set(USE_SEPVECTOR true)

if(BUILD_GCP)
  set(USE_GCP ON)
  set(GEN_IO_TYPES "${GEN_IO_TYPES} ,\"GCP_BUFFERS\"")
else(BUILD_GCP)

endif(BUILD_GCP)

if(SEPLIB_BUILD)

  set(GEN_IO_TYPES "${GEN_IO_TYPES} ,\"SEP\"")

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../sep3d
                      ${CMAKE_CURRENT_SOURCE_DIR}/../../sep
                      ${CMAKE_CURRENT_SOURCE_DIR}/../../include
                      ${CMAKE_CURRENT_BINARY_DIR}/../../)
  set(EXTERNAL_TARGETS ${EXTERNAL_TARGETS} sep3d sep)
  set(USE_SEP true)
else(SEPLIB_BUILD)
  if(BUILD_SEP)
    find_package(SEPlib)
  endif(BUILD_SEP)
  set(EXTERNAL_TARGETS)
  if(SEPlib_FOUND)
    message("INSTALLING WITH SEPLIB  DIR ${SEPlib_INCLUDE_DIRS}")
    include_directories(${SEPlib_INCLUDE_DIRS})
    set(EXTERNAL_TARGETS ${EXTERNAL_TARGETS} SEP:sep3d)
    set(USE_SEP true)
  else()
    message("NOT INSTALLING WITH SEPlib")
  endif()
endif(SEPLIB_BUILD)

if(USE_SEP)
  set(DEFAULTIO "SEP")
elseif(USE_RSF)
  set(DEFAULTIO "RSF")
else(USE_SEP)
  set(DEFAULTIO "JSON")
endif(USE_SEP)

configure_file("${PROJECT_SOURCE_DIR}/ioConfig.h.in"
               "${PROJECT_BINARY_DIR}/ioConfig.h")
set(genericIO_DIRNAME_lib "lib")
set(genericIO_DIRNAME_include "include")
set(genericIO_DIRNAME_cmake "cmake")
# cmake file
configure_file(genericIOConfig.cmake.in genericIOConfig.cmake @ONLY)

set(SUBPACKAGE true)

if(NOT TARGET sepVector)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/sepVector)
endif(NOT TARGET sepVector)

if(NOT TARGET buffers)
add_subdirectory(external/buffers)
endif(NOT TARGET buffers)


if(BUILD_TEST)
  enable_testing()
else()

endif(BUILD_TEST)

add_subdirectory(lib)

add_subdirectory(python)
if(BUILD_EXAMPLES)

add_subdirectory(examples)
endif(BUILD_EXAMPLES)

if(BUILD_UTIL)
  add_subdirectory(util)
endif(BUILD_UTIL)
  add_subdirectory(seis)

if(NOT SUPPORT_LIBRARIES_BUILD)
  install(FILES ${PROJECT_BINARY_DIR}/ioConfig.h DESTINATION include)
  add_subdirectory(env-modules)
endif(NOT SUPPORT_LIBRARIES_BUILD)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  genericIOConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/genericIOConfig.cmake"
  INSTALL_DESTINATION
  cmake/SEP/IO)

if(BUILD_DOC)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the
    # application
    add_custom_target(doc_doxygen_genericio ALL
                      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)
  else(DOXYGEN_FOUND)
    message(
      "Doxygen need to be installed to generate the doxygen documentation")
  endif(DOXYGEN_FOUND)
endif(BUILD_DOC)

# install(EXPORT genericIOTargets FILE genericIOTargets.cmake NAMESPACE SEP::IO
# DESTINATION cmake/sep/io)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${genericIO_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${genericIO_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_NAME "genericIO")
set(CPACK_PACKAGE_CONTACT "Robert Clapp")
include(CPack)
