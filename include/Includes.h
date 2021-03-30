#pragma once

// use STL library exceptions:
#include <exception>

// use Eigen library dense array/matrix objects:
#include <eigen/Eigen/Dense>

namespace Eigen {
  using Array23d  = Array<double, 2, 3>;
  using Array32d  = Array<double, 3, 2>;
  using Matrix23d = Matrix<double, 2, 3>;
  using Matrix32d = Matrix<double, 3, 2>;
  using Matrix13d = Matrix<double, 1, 3>;
  using Matrix31d = Matrix<double, 3, 1>;
}

// define possible boundary types
enum class boundaries : int { SOLID_WALL = -1, FREE_FLOW = -2, PERIODIC = -3, CUSTOM = -4 };

// define alias for idxing
using idx = Eigen::Index;
using Array = Eigen::Array3d;
using Storage1d = Eigen::ArrayXd;
using Storage3d = Eigen::Array3Xd;

// define computational precision
constexpr inline double tol = 1e-13;
