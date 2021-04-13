#pragma once
#include <Includes.h>
#include <Exceptions.h>

using Point = Array<2>;
using PointArray = Storage<2>;

//Simple operations on 2D Points
double len(Point const& a, Point const& b) noexcept;
double det(Point const& a, Point const& b) noexcept;
double triang_area(
  Point const& a,
  Point const& b,
  Point const& c) noexcept;
Point intersection(
  Point const& a1, Point const& a2,
  Point const& b1, Point const& b2);

double triang_average(Point const& p0, Point const& p1, Point const& p2,
  std::function<double(Point const&)> const& f);

template <class function>
double bisection(function const & f,
  double xmin = 0., double xmax = 1.,
  const double accuracy = 1e-14) {
  if (std::signbit(f(xmin)) == std::signbit(f(xmax)))
    return (abs(f(xmin)) < abs(f(xmax))) ? xmin : xmax;
  long N = static_cast<long>(log2(xmax - xmin) - log2(accuracy));
  double x;
  for (long i = 0; i <= N; i++) {
    x = 0.5*(xmin + xmax);
    ((std::signbit(f(xmin)) != std::signbit(f(x))) ? xmax : xmin) = x;
  }
  return x;
}

Eigen::Matrix32d gradient_coefs(const Eigen::Matrix32d& P);
