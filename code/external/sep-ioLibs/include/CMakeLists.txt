cmake_minimum_required(VERSION 2.8)

set(CMAKE_MACOSX_RPATH 1)

set(INCLUDES    
    sepcube.h	
   sep_old_external.h
	   sep_par.h	
prototypes.h  	   sep_pars_external.h	
	
sep3dc.h      sep_fortran.h	   		superset.h
sep3d.h       seplib.h		
      sep_main_external.h)

set(SCRIPTS sep_doc_it)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/su.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/su.h"
  )

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/segy.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/segy.h"
  )

install(FILES ${INCLUDES} DESTINATION include)
