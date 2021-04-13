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


