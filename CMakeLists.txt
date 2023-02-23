cmake_minimum_required(VERSION 3.15)
project(ksp)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_BUILD_TYPE_BACKUP ${CMAKE_BUILD_TYPE})
set(CMAKE_CXX_FLAGS_BACKUP ${CMAKE_CXX_FLAGS})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")

add_subdirectory(gtest)
add_subdirectory(gbench)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_BACKUP})
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_BACKUP})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
-fpermissive \
-Wall \
-Wextra \
-Wshadow=local \
-Wnon-virtual-dtor \
-pedantic \
-Wold-style-cast \
-Wuseless-cast \
-Wdouble-promotion \
-Wcast-align \
-Wunused \
-Woverloaded-virtual \
-Wconversion \
-Wsign-conversion \
-Wmisleading-indentation \
-Wduplicated-cond \
-Wduplicated-branches \
-Wlogical-op \
-Wnull-dereference \
-Wformat=2 \
-fno-rtti \
-march=native \
-pipe \
-fconcepts-diagnostics-depth=2 \
")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -fno-omit-frame-pointer -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
elseif(CMAKE_BUILD_TYPE MATCHES "Experiment")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -fno-exceptions") # todo move -fno-exceptions to default flags?
else(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g0 -DNDEBUG")
endif()

include(FindOpenMP)
if(NOT OPENMP_FOUND)
    message("OpenMP not found. Trying '-fopenmp=libiomp5'")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libiomp5")
else()
    message("OpenMP ${OpenMP_CXX_VERSION} found: ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(NOT DEFINED THREADS_KSP)
    set(THREADS_KSP 1)
endif()

if(NOT DEFINED THREADS_DS)
    set(THREADS_DS 1)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTHREADS_KSP=${THREADS_KSP} -DTHREADS_DS=${THREADS_DS}")

message("Compiler: ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}")
message("Flags: ${CMAKE_CXX_FLAGS}")

message("KSP Threads: ${THREADS_KSP}")
message("DS Threads: ${THREADS_DS}")

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/lib)

enable_testing()

add_subdirectory(tests)
add_subdirectory(benchmarks)
add_subdirectory(experiments)
add_subdirectory(tools)

