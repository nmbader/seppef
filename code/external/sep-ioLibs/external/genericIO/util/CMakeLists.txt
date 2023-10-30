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






  set(UTILSZ Compress  Decompress)
  foreach(util ${UTILSZ})
      ADD_EXECUTABLE(${util} ${util}.cc)
      TARGET_LINK_LIBRARIES(${util} genericCpp sepVector buffers hypercube tbb )
      install(TARGETS ${util} DESTINATION bin)
  endforeach(util)


install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/SEP2JSON DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/JSONIZE DESTINATION bin)




install(PROGRAMS ${progs} DESTINATION bin)

add_subdirectory(unix)
add_subdirectory(info)