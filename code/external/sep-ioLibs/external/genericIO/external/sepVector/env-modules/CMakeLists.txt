cmake_minimum_required(VERSION 2.8)
project(sepvector-env-modules)



  set(min  "${CMAKE_CURRENT_SOURCE_DIR}/sepVector.in")
  set(mout "${CMAKE_CURRENT_BINARY_DIR}/sepVector")
  configure_file (${min} ${mout})


install(FILES  ${mout}  DESTINATION share/modules)

