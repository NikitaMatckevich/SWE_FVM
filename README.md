SWE_FVM is a tiny finite-volume framework for solving 2D nonlinear shallow water equations with
source terms. It extensively exploits the ideas from the article [Liu et al.,
2018](https://www.sciencedirect.com/science/article/pii/S0021999118304996) to
build a robust positivity-preserving well-balanced interpolation of the unknown
fields.

To build the project, you will need:

* C/C++ compiler that supports the C++17 standard;
* Cmake version >= 3.8

Simply clone this repository to your machine and run Cmake:

$ mkdir build
$ cd build 
$ cmake -DCMAKE_BUILD_TYPE=Debug .. # or =Release
$ make

You will get the executable called SWE_FVM in your build directory. Change the
behavior of the program as you want by modyfying the src/Main.cpp.

Program supports the input from .ini files. Simply add something like
"config.ini" file into your project.

New flux discretizations, as well as time discretizations, can be added to the framework very simply (see the files
src/Fluxes.cpp and src/Solvers.cpp, respectively).

The author of the code apologizes for the lack of documentation. The project
will switch to Doxygen soon.
 
