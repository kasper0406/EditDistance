cmake_minimum_required (VERSION 2.6)
project (EditDistance)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/vendor/include/" "/usr/include/")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/" "/usr/lib/")

string(COMPARE EQUAL ${CMAKE_SYSTEM_NAME} "Linux" LINUX)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -flto -fasm-blocks -march=core-avx2 -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")

add_library(
  ipcm
  STATIC
  EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/client_bw.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/cpucounters.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/cpucounterstest.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/msr.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pci.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-memory.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-msr.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-numa.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-pcie.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-power.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-sensor.o
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/lib/pcm-tsx.o
)

set_target_properties(
  ipcm
  PROPERTIES
  LINKER_LANGUAGE CXX
)

set(ROOT_SRCS "main.cpp")
set(SIMPLE_SRCS "simple/simple.hpp")
set(COMPRESSION_SRCS "compression/SLP.hpp" "compression/DIST.hpp" "compression/slpalign.hpp" "compression/maxmultiply.hpp")
set(UTILS_SRCS "utils/allocator.hpp" "utils/matrix.hpp" "utils/math.hpp" "utils/unionfind.hpp")
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

target_link_libraries(EditDistance sdsl divsufsort divsufsort64 PcmMsr ipcm)
