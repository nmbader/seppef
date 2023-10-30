project(Generic-IO)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

#find_package(Boost 1.40 COMPONENTS program_options unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# This is a CMake example for Python


  pybind11_add_module(pyGenericIO MODULE pyIO.cpp)
  target_include_directories(
    pyGenericIO
    PUBLIC ${PYTHON_INCLUDE_DIR}
           ${pybind11_INCLUDE_DIR}
           ${PYTHON_INCLUDE_PATH}
           $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}../src>)

  target_link_libraries(pyGenericIO PRIVATE genericCpp)
  install(TARGETS pyGenericIO DESTINATION lib/python EXPORT genericIOTargets)
  install(FILES GenJob.py GenSplit.py genericIO.py DESTINATION lib/python)

set(PYSEARCH ${CMAKE_CURRENT_BINARY_DIR}/../lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../external/sepVector/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../external/sepVector/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../external/sepVector/lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../external/sepVector/external/hypercube/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../external/sepVector/external/hypercube/python)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_BINARY_DIR}/../external/sepVector/external/hypercube/lib)
set(PYSEARCH ${PYSEARCH}:${CMAKE_CURRENT_SOURCE_DIR}/../external/sepVector/external/python-solver/GenericSolver/python)

  if(BUILD_TEST)
  set(TESTS  Generic_IO GenSplit GenJob)
FOREACH(prg ${TESTS})
  add_test(NAME TEST_${prg}
  COMMAND ${CMAKE_COMMAND} -E env
  PYTHONPATH=${PYSEARCH}
  LD_LIBRARY_PATH=${PYSEARCH}
  DATAPATH=/tmp/
  ${PYTHON_EXECUTABLE} -m test_${prg}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDFOREACH(prg)
endif(BUILD_TEST)
