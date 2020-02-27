#pragma once
#include "ConfigParser.hpp" 

enum class scales { height, length, velocity, source, time };

class DimensionManager {
private:
  double h0, l0, c0;
  const double g = 9.8;
public:
  DimensionManager (Parser const & Parser) {
    h0 = Parser.get<double>("Common", "vertical_measure");
    l0 = Parser.get<double>("Common", "horizontal_measure");
    c0 = sqrt(g*h0);
  }
  template <scales T> double scale(double x) const {
    switch (T) {
      case scales::height   : return x * h0;
      case scales::length   : return x * l0;
      case scales::velocity : return x * c0;
      case scales::source   : return x * c0 / l0;
      case scales::time     : return x * l0 / c0;
    default: {
      std::ostringstream message;
      message << "Unexpected dimension in input";
      throw std::runtime_error(message.str());
    };
    }
  }
  template <scales T> double unscale(double x) const {
    switch (T) {
      case scales::height   : return x / h0;
      case scales::length   : return x / l0;
      case scales::velocity : return x / c0;
      case scales::source   : return x * l0 / c0;
      case scales::time     : return x * c0 / l0;
      default: {
        std::ostringstream message;
        message << "Unexpected dimension in input";
        throw std::runtime_error(message.str());
      };
    }
  }
};