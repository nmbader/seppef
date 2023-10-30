cmake_minimum_required(VERSION 3.5)
project(Generic-cpp)

#find_package(Boost 1.40 COMPONENTS program_options unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

set(SOURCES basicIO.cpp genericFile.cpp dictParams.cpp jsonGenericRegFile.cpp)
set(SOURCES
    ${SOURCES}
    jsonParamObj.cpp
    paramObj.cpp
    genericIO.cpp
    jsonGenericIrregFile.cpp
    buffersIrregFile.cpp
    memoryIO.cpp
    memoryFile.cpp)
set(SOURCES
    ${SOURCES}
    jsonGenericIO.cpp
    ioModes.cpp
    segyIO.cpp
    fileBuffersIO.cpp
    fileBuffersRegFile.cpp
    buffersRegFile.cpp)
set(INCLUDES genericIO.h jsonGenericRegFile.h jsonParamObj.h basicIO.h)
set(INCLUDES
    ${INCLUDES}
    jsonGenericIO.h
    paramObj.h
    genericFile.h
    ioModes.h
    segyIO.h
    dictParams.h
    fileBuffersIO.h
    fileBuffersRegFile.h
    memoryIO.h
    memoryFile.h
    buffersRegFile.h
    jsonGenericIrregFile.h
    buffersIrregFile.h)

if(USE_SEP)
  include_directories(${SEPlib_INCLUDE_DIRS})
  set(SOURCES
      ${SOURCES}
      sepIO.cpp
      sep3dFile.cpp
      sep_reg_file.cpp
      sep_params.cpp)
  set(INCLUDES ${INCLUDES} sep_params.h sepIO.h sep3dFile.h sep_reg_file.h)
endif(USE_SEP)

if(USE_RSF)
  include_directories(${RSFINCDIR})
  set(SOURCES ${SOURCES} rsfIO.cpp rsfParams.cpp rsfRegFile.cpp)
  set(INCLUDES ${INCLUDES} rsfIO.h rsfParams.h rsfRegFiles.h)
endif(USE_RSF)

if(BUILD_GCP)
  set(SOURCES ${SOURCES} gcpBuffersIO.cpp gcpBuffersRegFile.cpp)
  set(INCLUDES ${INCLUDE} gcpBuffersIO.h gcpBuffersRegFile.h)
endif(BUILD_GCP)

add_library(genericCpp ${SOURCES} ${INCLUDES})

target_include_directories(
  genericCpp
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(genericCpp sepVector buffers ${EXTERNAL_TARGETS})

if(NOT SUPPORT_LIBRARIES_BUILD)
  install(FILES ${INCLUDES} DESTINATION include)
  install(TARGETS genericCpp DESTINATION lib EXPORT genericIOExports)
endif(NOT SUPPORT_LIBRARIES_BUILD)

if(BUILD_TEST)
  add_executable(testMemoryFile test_memoryFile.cpp)
  target_link_libraries(testMemoryFile genericCpp  gtest gtest_main gtest )
  add_test(NAME memory_test COMMAND testMemoryFile)

  add_executable(testDict test_dictParams.cpp)
  target_link_libraries(testDict genericCpp  gtest gtest_main gtest )
  add_test(NAME dict_params COMMAND testDict)

  if(BUILD_GCP)
    add_executable(testGCP testGCP.cpp)
    target_link_libraries(testGCP genericCpp  gtest gtest_main gtest )
    add_test(NAME testGCP COMMAND testGCP)
  endif(BUILD_GCP)
endif(BUILD_TEST)

if(USE_SEP)

endif(USE_SEP)

if(BUILD_GCP)

endif(BUILD_GCP)
