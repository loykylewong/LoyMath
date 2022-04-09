cmake_minimum_required(VERSION 3.0.0)
project(loymath VERSION 0.1.0)

include(CTest)
enable_testing()

file(GLOB_RECURSE SRCS main.cpp test/*.cpp)
add_executable(loymath ${SRCS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
include(CPack)
