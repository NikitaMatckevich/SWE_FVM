#pragma once

// use Eigen library dense array/matrix objects:
#include "eigen/Eigen/Dense"

namespace Eigen {
  using Array23d  = Array<double, 2, 3>;
  using Array32d  = Array<double, 3, 2>;
  using Matrix23d = Matrix<double, 2, 3>;
  using Matrix32d = Matrix<double, 3, 2>;
  using Matrix13d = Matrix<double, 1, 3>;
  using Matrix31d = Matrix<double, 3, 1>;
  using Matrix12d = Matrix<double, 1, 2>;
  using Matrix21d = Matrix<double, 2, 1>;
}

// define possible boundary types
enum class Boundaries : int { SOLID_WALL = -1, FREE_FLOW = -2, PERIODIC = -3, CUSTOM = -4 };

// define alias for idxing
using Idx = Eigen::Index;

template <size_t k>
using Array = Eigen::Array<double, k, 1>;

template <size_t k>
using Storage = Eigen::Array<double, k, Eigen::Dynamic>;

// define computational precision
constexpr inline double tol = 1e-13;
