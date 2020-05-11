cmake_minimum_required (VERSION 3.8)

project ("SWE_FVM")
set(CMAKE_CXX_STANDARD 17)

add_executable( SWE_FVM "TestAll.cpp"
						"ConfigParser.cpp"
						"DimensionManager.cpp"
						"PointOperations.cpp"
						"TriangMesh.cpp"
						"Bathymetry.cpp"
						"ConsAssigner.cpp"
						"ValueField.cpp"
						"MUSCLObject.cpp"
						"Solver.cpp"
						"KurganovFluxMethod.cpp")