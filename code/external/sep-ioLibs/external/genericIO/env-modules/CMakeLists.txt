cmake_minimum_required(VERSION 2.8)
project(genericIO-envModules)



  set(min  "${CMAKE_CURRENT_SOURCE_DIR}/genericIO.in")
  set(mout "${CMAKE_CURRENT_BINARY_DIR}/genericIO")
  configure_file (${min} ${mout})


install(FILES  ${mout}  DESTINATION share/modules)

