#pragma once
#include "ConfigParser.h"
#include <Exceptions.h>

enum class Scales { height, length, velocity, source, time };

struct DimensionManager {

  explicit DimensionManager(const Parser& parser);
  
  template <Scales T>
  double Scale(double x) const {
    switch (T) {
    case Scales::height: return x * m_h0;
    case Scales::length: return x * m_l0;
    case Scales::velocity: return x * m_c0;
    case Scales::source: return x * m_c0 / m_l0;
    case Scales::time: return x * m_l0 / m_c0;
    default: throw DimensionError("invalid dimension type in scale function");
    }
  }

  template <Scales T>
  double Unscale(double x) const {
    switch (T) {
    case Scales::height: return x / m_h0;
    case Scales::length: return x / m_l0;
    case Scales::velocity: return x / m_c0;
    case Scales::source: return x * m_l0 / m_c0;
    case Scales::time: return x * m_c0 / m_l0;
    default: throw DimensionError("invalid dimension type in unscale function");
    }
  }

 private:

  double m_h0, m_l0, m_c0;
};

