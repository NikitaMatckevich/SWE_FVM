#include "PointOperations.h"

double length(Point const& a, Point const& b) noexcept {
  return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]));
}
double det(Point const& a, Point const& b) noexcept {
  return a[0] * b[1] - b[0] * a[1];
}
double triangle_area(
  Point const& a,
  Point const& b,
  Point const& c) noexcept
{
  return 0.5 * abs(det(b - a, c - a));
}
Point intersection(
  Point const& a1, Point const& a2,
  Point const& b1, Point const& b2)
{
  double den = det(a2 - a1, b2 - b1);
  if (den == 0.)
    throw std::runtime_error("trying to find intersection of parallel vectors");
  return a1 + (a2 - a1) * det(b1 - a1, b2 - b1) / den;
}

double avg_on_triangle(Point const& p0, Point const& p1, Point const& p2,
  std::function<double(Point const&)> const& f) {
  double S = triangle_area(p0, p1, p2);
  constexpr unsigned N = 500;
  constexpr double   h = 1./N;
  Point di = h * (p1 - p0);
  Point dj = h * (p2 - p0);
  Point dt = 1. / 3. * (di + dj);
  double sum = 0.;
  Point pi = p0;
  for (unsigned i = 0; i < N; i++) {
    Point pt = pi + dt;
    for (unsigned j = 0; j < N - i - 1; j++) {
      sum += h * f(pt);
      sum += h * f(pt + dt);
      pt += dj;
    }
    sum += h * f(pt);
    pi += di;
  }
  return h * sum;
}