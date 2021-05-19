__SWE_FVM__ is a tiny finite-volume framework for solving 2D nonlinear shallow water equations with
source terms. It extensively exploits the ideas from the article [Liu et al.,
2018](https://www.sciencedirect.com/science/article/pii/S0021999118304996) to
build a robust positivity-preserving well-balanced interpolation of the unknown
fields.

To build the project, you will need:

* C/C++ compiler that supports the C++17 standard;
* Cmake version >= 3.8.

Simply clone this repository to your machine and run Cmake:

```
$ mkdir build
$ cd build 
$ cmake -DCMAKE_BUILD_TYPE=Debug .. # or =Release
$ make
```

You will get the executable called SWE_FVM in your _build_ directory. Change the
behavior of the program as you want by modyfying the src/Main.cpp.

Program supports the input from .ini files.

New flux discretizations, as well as time discretizations, can be added to the framework very simply.

Please check the [documentation](https://nikitamatckevich.github.io/SWE_FVM/)
page for more details.
