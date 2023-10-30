cmake_minimum_required(VERSION 2.8)
set(CMAKE_MACOSX_RPATH 1)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep3d)


set (SOURCES 
axes.c	   drn.c	   header.c   list.c	    rite.c	   tagstream.c
conform.c  extra_keys.c    headval.c   rite_file.c    transfer.c
convert.c  fortran_junk.c  init.c     ntraces.c     rite.inc	   type_usage.c
coord.c    	    reed.c	   
debug.c    grid.c	   keys.c     struct_sepc.c)

if(CMAKE_Fortran_COMPILER)
set(SOURCES ${SOURCES} struct_sepf.f90)
endif(CMAKE_Fortran_COMPILER)



add_library(superset ${SOURCES})

target_include_directories(superset
	PRIVATE          
	   ${CMAKE_CURRENT_SOURCE_DIR}/../include 
	PUBLIC
	   $<INSTALL_INTERFACE:include> 
	   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
	   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>   
      )
install(TARGETS superset DESTINATION lib EXPORT  SEPlibTargets)
target_link_libraries(superset sep3d sep)
SET_TARGET_PROPERTIES(superset PROPERTIES LINKER_LANGUAGE Fortran)

if(SU_LIBRARIES)
include_directories(${SU_INCLUDES})
set(SOURCE2 reed_su.c rite_su.c init_su.c)
add_library(sepsu ${SOURCE2})
install(TARGETS sepsu DESTINATION lib)
endif(SU_LIBRARIES)
