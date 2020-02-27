#pragma once
#include <Eigen/Dense>
namespace Eigen {
  using Array23d = Array <double, 2, 3>;
  using Array32d = Array <double, 3, 2>;
  using Matrix23d = Matrix<double, 2, 3>;
  using Matrix32d = Matrix<double, 3, 2>;
}
using point = Eigen::Array2d;

//Simple operations on 2D points
inline double length(point const & a, point const & b) {
  return sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]));
}
inline double det(point const & a, point const & b) {
  return a[0]*b[1] - b[0]*a[1];
}
inline double triangle_area(point const & a,
                            point const & b,
                            point const & c) {
  return 0.5 * abs(det(b-a, c-a));
}
point intersection(point const & a1, point const & a2,
                   point const & b1, point const & b2) {
  return a1 + (a2 - a1) * det(b1 - a1, b2 - b1) / det(a2 - a1, b2 - b1);
}

template <typename T>
inline Eigen::Matrix<T, 3, 1> crcShiftL(Eigen::Matrix<T, 3, 1> const & v) {
  return { v[1], v[2], v[0] };
}
template <typename T>
inline Eigen::Matrix<T, 3, 1> crcShiftR(Eigen::Matrix<T, 3, 1> const & v) {
  return { v[2], v[0], v[1] };
}