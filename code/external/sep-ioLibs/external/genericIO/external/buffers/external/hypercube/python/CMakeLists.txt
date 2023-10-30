project(hypercube-python)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

find_package (Python3 COMPONENTS Interpreter  REQUIRED)



  pybind11_add_module(pySEPException MODULE pySEPException.cpp)
  target_link_libraries(pySEPException PRIVATE hypercube)
  target_include_directories(
    pySEPException
    PUBLIC ${PYTHON_INCLUDE_DIR}
           ${pybind11_INCLUDE_DIR}
           $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}../src>)

  pybind11_add_module(pyHypercube MODULE pyHyper.cpp)
  target_link_libraries(pyHypercube PRIVATE hypercube)
  target_include_directories(
    pyHypercube
    PUBLIC ${PYTHON_INCLUDE_DIR}
           ${pybind11_INCLUDE_DIR}
           $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}../src>)

  pybind11_add_module(pyDebug MODULE pyDebug.cpp)
  target_link_libraries(pyDebug PRIVATE hypercube)
  target_include_directories(
    pyDebug
    PUBLIC ${PYTHON_INCLUDE_DIR}
           ${pybind11_INCLUDE_DIR}
           $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}../src>)
  install(TARGETS pyHypercube pySEPException pyDebug
          DESTINATION lib/python
          EXPORT HypercubeTargets)
  install(FILES Hypercube.py DESTINATION lib/python)
if(Python3_Interpreter_FOUND)
  add_test(NAME hypercubePython COMMAND Python3::Interpreter -m pytest -r a -v  ${CMAKE_CURRENT_SOURCE_DIR})
endif()
