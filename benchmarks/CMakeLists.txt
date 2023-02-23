set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(GLOB BENCHMARKS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

foreach( sourcefile ${BENCHMARKS} )
    string( REPLACE ".cpp" "" name ${sourcefile} )
    add_executable( "${name}Benchmark" ${sourcefile} )
    target_link_libraries("${name}Benchmark" benchmark::benchmark)
endforeach()