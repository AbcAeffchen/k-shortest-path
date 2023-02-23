# Add an executable containing all tests

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

file(GLOB test_src "*.cpp")
add_executable(AllTests ${test_src})
target_link_libraries(AllTests gtest_main)