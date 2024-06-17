#pragma once
#include <stdexcept>

struct DomainError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

