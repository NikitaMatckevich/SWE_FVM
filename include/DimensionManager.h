#pragma once
#include <ConfigParser.h>
#include <Exceptions.h>

enum class scales { height, length, velocity, source, time };

class DimensionManager {
  Parser const& Parser_;
  double h0_, l0_, c0_;
public:
  DimensionManager(Parser const&);
  template <scales T>
  double scale(double x) const {
    switch (T) {
    case scales::height: return x * h0_;
    case scales::length: return x * l0_;
    case scales::velocity: return x * c0_;
    case scales::source: return x * c0_ / l0_;
    case scales::time: return x * l0_ / c0_;
    default: throw DimensionError("invalid dimension type in scale function");
    }
  }
  template <scales T>
  double unscale(double x) const {
    switch (T) {
    case scales::height: return x / h0_;
    case scales::length: return x / l0_;
    case scales::velocity: return x / c0_;
    case scales::source: return x * l0_ / c0_;
    case scales::time: return x * c0_ / l0_;
    default: throw DimensionError("invalid dimension type in unscale function");
    }
  }
};
