#pragma once
#include "Includes.h"

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct Solver {
  static inline constexpr double tol = 1e-13;
};