#pragma once

// STL library:
#include <exception>
#include <functional>
#include <algorithm>

// Eigen library:
#include "Eigen/Dense"
namespace Eigen {
  using Array23d = Array<double, 2, 3>;
  using Array32d = Array<double, 3, 2>;
  using Matrix23d = Matrix<double, 2, 3>;
  using Matrix32d = Matrix<double, 3, 2>;
  using Matrix13d = Matrix<double, 1, 3>;
  using Matrix31d = Matrix<double, 3, 1>;
}

enum boundaries { SOLID_WALL = -1, FREE_FLOW = -2, PERIODIC = -3, CUSTOM = -4 };

using index = Eigen::Index;

constexpr inline double tol = 1e-13;