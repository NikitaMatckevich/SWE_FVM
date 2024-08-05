cmake_minimum_required(VERSION 3.24)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Needed for vim code completion plugin
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fomit-frame-pointer")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Current build type: ${CMAKE_BUILD_TYPE}")

include(FetchContent)
FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(Eigen3)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

file(GLOB SOURCE_FILES "src/*[a-zA-Z].cpp")
file(GLOB PYTHON_FILES "pybind/*[a-zA-Z].cpp")
include_directories("include")

pybind11_add_module(SWE_FVM SHARED ${SOURCE_FILES} ${PYTHON_FILES})
target_link_libraries(SWE_FVM PUBLIC Eigen3::Eigen)
install(TARGETS SWE_FVM LIBRARY DESTINATION .)

add_executable(test examples/main.cpp)
target_link_libraries(test PRIVATE SWE_FVM)
message(STATUS "Test executable generated in ${RUNTIME_OUTPUT_DIRECTORY}")
