#Required CMake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

#We need to add to how to find TBB which sits in the cmake dir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


#Name of our cmake project
project(PEF LANGUAGES CXX)

#remove CMAKE_INSTALL_PREFIX. Use DESTDIR instead when installing software to a given directory
#example: make install DESTDIR="/usr/local/"
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "" CACHE PATH "Removing CMAKE_INSTALL_PREFIX. Use DESTDIR instead." FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


#Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Ofast")

find_package(FFTW REQUIRED)


add_subdirectory(external/sep-ioLibs)


#Include SIGPACK directory
#include_directories(/net/server2/homes/sep/nmbader/external/sigpack/sigpack)

#Include ARMADILLO directory
#include_directories(/net/server2/homes/sep/nmbader/external/armadillo/include)

#Adding path to link ARMADILLO library
#link_directories(/net/server2/homes/sep/nmbader/external/armadillo/lib64)

# We need Python for the figures
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_PATH})

#TBB and OpenMP required
find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#we also need boost
FIND_PACKAGE( Boost 1.40 COMPONENTS program_options unit_test_framework REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

#add the rules for building the C++ code (execute the CMakelists file in the subdirectory)
add_subdirectory(src)
