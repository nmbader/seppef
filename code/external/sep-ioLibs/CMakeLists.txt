if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif(POLICY CMP0048)
project(SEPlib-io-libs VERSION 9.0.0 LANGUAGES CXX C Fortran)

if(NOT SUBPACKAGE)
  option(BUILD_PYTHON "Build python interfaces" ON)
  option(BUILD_TEST "Build tests" OFF)
endif()

cmake_minimum_required(VERSION 3.2)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_POSITION_INDEPENDENT_CODE True)
#
# SEPLIB can be built with or without shared libraries. If shared libraries are
# desired, we suggest building with shared libraries first using
#
# cmake -DBUILD_SHARED_LIBS=True
#
# and then rebuilding without shared libraries using
#
# cmake -DBUILD_SHARED_LIBS=False
#
# so that the standard program executables are statically linked.
#
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(SEPlib-io_MAJOR_VERSION 9)
set(SEPlib-io_MINOR_VERSION 0)
set(SEPlib-io_VERSION_STRING
    ${SEPlib-io_MAJOR_VERSION}.${SEPlib-io_MINOR_VERSION})

set(SEPlib_VERSION ${SEPlib-io_VERSION_STRING})
enable_language(C OPTIONAL Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")

# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# Offer the user the choice of overriding the installation directories
# set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
# set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH

# SET(SEPlib_DIRNAME_lib     ${CMAKE_INSTALL_PREFIX}/lib)
# SET(SEPlib_DIRNAME_include ${CMAKE_INSTALL_PREFIX}/include)
# SET(SEPlib_DIRNAME_cmake ${CMAKE_INSTALL_PREFIX}/cmake)

if(NOT SUBPACKAGE)
  option(BUILD_TEST "Build with test support" OFF)
  option(BUILD_DEBUG_LIBS "Build gperftool libs" OFF)
  option(BUILD_GCP "Build GCP" OFF)

else(NOT SUBPACKAGE)

endif(NOT SUBPACKAGE)

find_package(SU)
if(SU_LIBRARIES)
  set(SU_SUPPORT true)
endif(SU_LIBRARIES)

include(CheckIncludeFiles)
check_include_files(termios.h HAVE_TERMIOS_H)
check_include_files(termio.h HAVE_TERMIO_H)
check_include_files(stdio.h HAVE_STDIO_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(sys/ioctl.h HAVE_SYS_IOCTL_H)
check_include_files(sys/file.h HAVE_SYS_FILE_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/tape.h HAVE_SYS_TAPE_H)
check_include_files(sys/uio.h HAVE_SYS_UIO_H)
check_include_files(sgtty.h HAVE_SGTTY_H)
check_include_files(signal.h HAVE_SIGNAL_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(sys/termios.h HAVE_SYS_TERMIOS_H)
check_include_files(sys/utsname.h HAVE_SYS_UTSNAME_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(fcntl.h HAVE_FCNTL_H)
check_include_files(dirent.h HAVE_DIRENT_H)
check_include_files(sys/ndir.h HAVE_SYS_NDIR_H)
check_include_files(sys/dir.h HAVE_SYS_DIR_H)
check_include_files(ndir.h HAVE_NDIR_H)
check_include_files(limits.h HAVE_LIMITS_H)
check_include_files(errno.h HAVE_ERRNO_H)
check_include_files(sys/socket.h HAVE_SYS_SOCKET_H)
check_include_files(rpc/rpc.h HAVE_RPC_RPC_H)
check_include_files(rpc/types.h HAVE_RPC_TYPES_H)
check_include_files(ctime.h HAVE_CTIME_H)
check_include_files(utime.h HAVE_UTIME_H)


configure_file("${PROJECT_SOURCE_DIR}/sepConfig.h.in"
               "${PROJECT_BINARY_DIR}/sepConfig.h")

#
set(SEPCONFIGDIR ${CMAKE_CURRENT_BINARY_DIR})
set(SEPLIB_BUILD yes)
if(NOT TARGET genericCpp AND NOT SUBPACKAGE)
    execute_process(COMMAND git submodule update --init --recursive --
                            external/genericIO
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif(NOT TARGET genericCpp AND NOT SUBPACKAGE)

set(SUBPACKAGE true)
add_subdirectory(sep)
add_subdirectory(sep3d)
add_subdirectory(superset)
add_subdirectory(include)

if(NOT TARGET genericCpp)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/genericIO)
endif(NOT TARGET genericCpp)


export(TARGETS sep sep3d superset
       FILE "${PROJECT_BINARY_DIR}/SEPlibTargets.cmake")
export(PACKAGE SEPlib)

include(CMakePackageConfigHelpers)
set(SEPlib_INCLUDE_DIR "include")

configure_file(SEPlibConfig.cmake.in
               "${CMAKE_CURRENT_BINARY_DIR}/SEPlibConfig.cmake" @ONLY)

# cmake file configure_file (SEPlibConfig.cmake.in SEPlibConfig.cmake @ONLY)
# configure_file (SEPlibConfigVersion.cmake.in SEPlibConfigVersion.cmake @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SEPlibConfig.cmake
        DESTINATION cmake/SEP/)

# install(FILES "${PROJECT_BINARY_DIR}/cmake/SEPlibConfig.cmake"
# "${PROJECT_BINARY_DIR}/cmake/SEPlibConfigVersion.cmake" DESTINATION
# "${SEPlib_DIRNAME_cmake}" COMPONENT dev)

install(EXPORT SEPlibTargets
        FILE SEPlibTargets.cmake
        NAMESPACE SEP:
        DESTINATION cmake/SEP)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sepConfig.h DESTINATION include)



include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${SEPlib-io_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${SEPlib-io_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_NAME "SEPlib-io")
set(CPACK_PACKAGE_CONTACT "Robert Clapp")
include(CPack)
