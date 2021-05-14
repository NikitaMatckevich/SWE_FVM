SWE_FVM is a tiny finite-volume framework for solving 2D nonlinear shallow water equations with
source terms. It extensively exploits the ideas from the article [Liu et al.,
2018](https://www.sciencedirect.com/science/article/pii/S0021999118304996) to
build a robust positivity-preserving well-balanced interpolation of the unknown
fields. Flux discretizations, as well as time discretizations, can be added to the framework very simply (see the files
src/Fluxes.cpp and src/Solvers.cpp, respectively).

The author of the code apologizes for the lack of documentation. The project
will switch to Doxygen soon.
 
