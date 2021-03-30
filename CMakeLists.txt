cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Needed for vim code completion plugin

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pthread -DDEBUG -g -Wall -Wextra -pedantic -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -pthread -fomit-frame-pointer")

message(STATUS "Current build type: ${CMAKE_BUILD_TYPE}")

project ("SWE_FVM")

#file(GLOB SRC "src/*[a-zA-Z].cpp")
set(SRC "src")
include_directories("include")

add_executable(SWE_FVM
	${SRC}/Main.cpp;
	${SRC}/ConfigParser.cpp;
	${SRC}/DimensionManager.cpp;
	${SRC}/PointOperations.cpp;
	${SRC}/TriangMesh.cpp;
  ${SRC}/Bathymetry.cpp
	${SRC}/ConsAssigner.cpp;
	${SRC}/ValueField.cpp;
  ${SRC}/MUSCLObject.cpp;
  ${SRC}/Solver.cpp;
	${SRC}/SpaceIntegrators.cpp)
