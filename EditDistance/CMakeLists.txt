cmake_minimum_required (VERSION 2.6)
project (EditDistance)

string(COMPARE EQUAL ${CMAKE_SYSTEM_NAME} "Linux" LINUX)

IF(APPLE)
  set(CMAKE_CXX_COMPILER "clang++")
#  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
ENDIF(APPLE)

set(ROOT_SRCS "main.cpp")
set(SIMPLE_SRCS "simple/simple.hpp")
set(COMPRESSION_SRCS "compression/SLP.hpp" "compression/DIST.hpp" "compression/slpalign.hpp" "compression/maxmultiply.hpp")
set(UTILS_SRCS "utils/allocator.hpp" "utils/matrix.hpp" "utils/math.hpp")
set(BENCHMARK_SRCS "benchmark/benchmarker.hpp")
set(TEST_SRCS "test/tests.hpp")

option(USE_COUNTING_ALLOCATOR "Use allocator counting memory consumption." ON)
if(USE_COUNTING_ALLOCATOR)
  add_definitions(-DUSE_COUNTING_ALLOCATOR)
endif()

add_executable(EditDistance ${ROOT_SRCS}
                            ${UTILS_SRCS}
                            ${SIMPLE_SRCS}
                            ${COMPRESSION_SRCS}
                            ${BENCHMARK_SRCS}
                            ${TEST_SRCS})

source_group("Simple"      FILES ${SIMPLE_SRCS})
source_group("Compression" FILES ${COMPRESSION_SRCS})
source_group("Utils"       FILES ${UTILS_SRCS})
source_group("Benchmark"   FILES ${BENCHMARK_SRCS})
source_group("Tests"       FILES ${TEST_SRCS})
