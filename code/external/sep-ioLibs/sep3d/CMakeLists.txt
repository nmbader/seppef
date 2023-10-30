cmake_minimum_required(VERSION 3.2)
set(CMAKE_MACOSX_RPATH 1)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep)

set(SOURCES
    copy_grid.c
    getkeyfirstbyte.c
    getvalheaders.c
    putnumberkeys.c
    copykeys.c
    getkeyfmt.c
    handle_file_pointers.c
    putvalbindex.c
    extractvalbyindex.c
    getkeyindex.c
    init_3d.c
    putvalbname.c
    fortran_junk.c
    getkeyname.c
    insertvalbyindex.c
    putvalheaders.c
    getkeytype.c
    putaxis.c
    reorder2.c
    getaxis.c
    getnumberaxes.c
    putgetgridwindow.c
    reorder.c
    getgridtag.c
    getnumberkeys.c
    putgridtag.c
    sep_3d_close.c
    getheaderbytes.c
    getvalbindex.c
    putheadertag.c
    tag_is_pipe.c
    getheadertag.c
    getvalbname.c
    putkey.c)

if(CMAKE_Fortran_COMPILER)
  set(SOURCES ${SOURCES} sep3d.f90)
endif(CMAKE_Fortran_COMPILER)
add_library(sep3d ${SOURCES})
target_link_libraries(sep3d sep)

SET_TARGET_PROPERTIES(sep3d PROPERTIES LINKER_LANGUAGE Fortran)
target_include_directories(
  sep3d
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(TARGETS sep3d DESTINATION lib EXPORT SEPlibTargets)
