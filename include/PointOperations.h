#pragma once
#include <Exceptions.h>
#include <Includes.h>

using Point = Array<2>;
using PointArray = Storage<2>;

//Simple operations on 2D Points
double Len(const Point& a, const Point& b) noexcept;
double Det(const Point& a, const Point& b) noexcept;

double TriangArea(
  const Point& a,
  const Point& b,
  const Point& c) noexcept;

Point Intersection(
  const Point& a1, const Point& a2,
  const Point& b1, const Point& b2);

double TriangAverage(const Point& p0, const Point& p1, const Point& p2,
  const std::function<double(const Point&)>& f);

template <class Function>
double Bisection(
  const Function& f,
  double xmin = 0., double xmax = 1.,
  const double accuracy = 1e-14)
{
  if (std::signbit(f(xmin)) == std::signbit(f(xmax)))
    return (abs(f(xmin)) < abs(f(xmax))) ? xmin : xmax;
  long n = static_cast<long>(log2(xmax - xmin) - log2(accuracy));
  double x;
  for (long i = 0; i <= n; i++) {
    x = 0.5*(xmin + xmax);
    ((std::signbit(f(xmin)) != std::signbit(f(x))) ? xmax : xmin) = x;
  }
  return x;
}

Eigen::Matrix32d GradientCoefs(const Eigen::Matrix32d& p);
