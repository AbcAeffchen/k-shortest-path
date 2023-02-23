set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(GLOB_RECURSE EXPERIMENTS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/code" *.cpp)

foreach( sourcefile ${EXPERIMENTS} )
    message(${sourcefile})
    string( REPLACE ".cpp" "" name ${sourcefile} )
    add_executable( ${name} code/${sourcefile} )
endforeach()