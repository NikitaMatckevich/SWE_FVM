#include "Solver.h"

using namespace utils;

template <class FluxMethod>
struct EulerForwardSolver : BaseSolver<FluxMethod> {
  using BaseSolver::BaseSolver;
  void step(double dt) {
    reconstruct_all();
    const auto& m = bathymetry().mesh();
    for (size_t i = 0; i < m.num_triangles(); ++i)
      Vol_.cons(i) += dt * rhs(i);
  }
};