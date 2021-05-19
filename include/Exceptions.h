#pragma once
#include <stdexcept>

struct ParserError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct DimensionError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct MeshError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

