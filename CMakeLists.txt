cmake_minimum_required (VERSION 3.8)

project ("SWE_FVM")
set(CMAKE_CXX_STANDARD 17)

add_executable( SWE_FVM "Main.cpp"
						"ConfigParser.cpp"
						"DimensionManager.cpp"
						"PointOperations.cpp"
						"TriangMesh.cpp"
						"Bathymetry.cpp"
						"ValueFieldAssign.cpp"
						"ValueFieldUpdate.cpp")