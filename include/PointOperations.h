#pragma once
#include <Exceptions.h>
#include <Includes.h>

using Point = Array<3>;
using PointArray = Storage<3>;

double Len(const Point& a, const Point& b) noexcept;
double Det(const Point& a, const Point& b) noexcept;

double TriangArea(
  const Point& a,
  const Point& b,
  const Point& c) noexcept;

Point Intersection(
  const Point& a1, const Point& a2,
  const Point& b1, const Point& b2);

template <size_t k, size_t n>
Array<k> TriangAverage(const Point& p0, const Point& p1, const Point& p2,
  const std::function<Array<k>(const Point&)>& f)
{
  constexpr double h = 1./n;
  const Point di = h * (p1 - p0);
  const Point dj = h * (p2 - p0);
  const Point dt = 1. / 3. * (di + dj);
  Array<k> sum = Array<k>::Zero();
  
  Point pi = p0;
  
  for (unsigned i = 0; i < n; i++) {
    Point pt = pi + dt;
    for (unsigned j = 0; j < n - i - 1; j++) {
      sum += h * f(pt);
      sum += h * f(pt + dt);
      pt += dj;
    }
    sum += h * f(pt);
    pi += di;
  }

  return h * sum;
}

double Bisection(
  const std::function<double(double)>& f,
  double xmin = 0., double xmax = 1.,
  const int accuracy = 50);

Eigen::Vector2d Gradient(const Eigen::Matrix3d& p);
