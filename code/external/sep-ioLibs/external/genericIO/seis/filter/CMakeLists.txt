project(Generic-IO-seis-model)
cmake_minimum_required(VERSION 2.8)





set(PYSEARCH ${CMAKE_CURRENT_BINARY_DIR}/../../lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../../python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../../external/sepVector/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../../external/sepVector/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../../external/sepVector/lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../../external/sepVector/external/hypercube/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../../external/sepVector/external/hypercube/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../../external/sepVector/external/hypercube/lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../../external/sepVector/external/python-solver/GenericSolver/python)

  if(BUILD_TEST)
  set(TESTS   )
FOREACH(prg ${TESTS})
  add_test(NAME TEST_${prog}
  COMMAND ${CMAKE_COMMAND} -E env
  PYTHONPATH=${PYSEARCH}
  LD_LIBRARY_PATH=${PYSEARCH}
  DATAPATH=/tmp/
  ${PYTHON_EXECUTABLE} -m test_${prg}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDFOREACH(prg)
endif(BUILD_TEST)

set(SCRIPTS  GenInterp)
FOREACH(prg ${SCRIPTS})
   configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/${prg}.py"
  "${CMAKE_CURRENT_BINARY_DIR}/${prg}"
  @ONLY
  )
   install(PROGRAMS  ${CMAKE_CURRENT_BINARY_DIR}/${prg} DESTINATION bin)
ENDFOREACH(prg)

install(PROGRAMS ${progs} DESTINATION bin)