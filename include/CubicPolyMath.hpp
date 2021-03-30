#pragma once
#include <iostream>
#include <cassert>
#include <functional>

namespace cubic {
  class poly {
  private:
    double b, c, d;
  public:
    poly(double p, double q) // x^3 + p*x + q
      : b(0)
      , c(p)
      , d(q) {}
    poly(double b, double c, double d) // x^3 + b*x^2 + c*x + d
      : b(b)
      , c(c)
      , d(d) {}
    inline double operator()(double x) const {
      return x*x*x + b*x*x + c*x + d;
    }
  };
}

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
