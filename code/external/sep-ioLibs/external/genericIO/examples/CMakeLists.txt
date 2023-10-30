project(Generic-IO-util)
cmake_minimum_required(VERSION 2.8)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cpp)

FIND_PACKAGE(TBB REQUIRED)


# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

FOREACH(lib ${GENERIC_EXPORT_LIB_DIRS})
  link_directories(${lib})
ENDFOREACH(lib)
include(ExternalProject)


  externalproject_add(
    gperftools_src
    PREFIX
    gperftools_src
    INSTALL_DIR
    ${CMAKE_INSTALL_PREFIX}
    URL
    https://github.com/gperftools/gperftools/releases/download/gperftools-2.7/gperftools-2.7.tar.gz
    CONFIGURE_COMMAND
    ./configure
    --enable-frame-pointers
    --prefix=<INSTALL_DIR>
    BUILD_IN_SOURCE
    1)

  externalproject_get_property(gperftools_src INSTALL_DIR)
  set(GPERF_LIBDIR ${INSTALL_DIR}/lib)
  set(PROFILER_LIB profiler)
  set(TCMALLOC_LIB tcmalloc)
  add_library(${PROFILER_LIB} UNKNOWN IMPORTED)
  add_library(${TCMALLOC_LIB} UNKNOWN IMPORTED)
  set_property(
    TARGET ${PROFILER_LIB}
    PROPERTY
      IMPORTED_LOCATION
      ${GPERF_LIBDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}profiler${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
  set_property(
    TARGET ${TCMALLOC_LIB}
    PROPERTY
      IMPORTED_LOCATION
      ${GPERF_LIBDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}tcmalloc${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
  add_dependencies(${PROFILER_LIB} gperftools_src)
  add_dependencies(${TCMALLOC_LIB} gperftools_src)



if(BUILD_GCP)
set(EXAMPlES ${EXAMPLES} ReadWriteGCP)
endif()

foreach(exam ${EXAMPlES})
  ADD_EXECUTABLE(${exam} ${exam}.cc)
  TARGET_LINK_LIBRARIES(${exam} genericCpp hypercube)
  foreach(lib ${GENERIC_EXPORT_LIBS})
    TARGET_LINK_LIBRARIES(${exam} ${lib})
  endforeach(lib)
  TARGET_LINK_LIBRARIES(${exam} sepVector buffers tbb ${PROFILER_LIB} ${TCMALLOC_LIB})
  install(TARGETS ${exam} DESTINATION bin)
endforeach(exam)

#find_package(ZFP)





