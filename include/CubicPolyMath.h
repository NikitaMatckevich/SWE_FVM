#pragma once
#include <cassert>
#include <functional>
#include <iostream>

struct CubicPoly {

  CubicPoly(double d = 0, double c = 0, double b = 0) // x^3 + b*x^2 + c*x + d
    : m_b(b)
    , m_c(c)
    , m_d(d) {}

  inline double operator()(double x) const {
    return x*x*x + m_b*x*x + m_c*x + m_d;
  }

 private:
  double m_b, m_c, m_d;
};
