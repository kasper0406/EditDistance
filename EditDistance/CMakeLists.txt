cmake_minimum_required (VERSION 2.6)
project (EditDistance)

string(COMPARE EQUAL ${CMAKE_SYSTEM_NAME} "Linux" LINUX)

IF(APPLE)
  set(CMAKE_CXX_COMPILER "clang++")
#  set(CMAKE_BUILD_TYPE Release)
#  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
ENDIF(APPLE)

add_executable(EditDistance main.cpp)