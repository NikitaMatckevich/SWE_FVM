cmake_minimum_required(VERSION 3.24)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Needed for vim code completion plugin
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g -Wall -Wextra -pedantic -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fomit-frame-pointer")

FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(Eigen3)

message(STATUS "Current build type: ${CMAKE_BUILD_TYPE}")

project("SWE_FVM")

file(GLOB SRC "src/*[a-zA-Z].cpp")
message("Sources:")
message(${SRC})

include_directories("include")

add_executable(SWE_FVM ${SRC})
target_link_libraries(SWE_FVM PRIVATE Eigen3::Eigen) 
