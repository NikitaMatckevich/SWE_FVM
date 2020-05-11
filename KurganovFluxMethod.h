#pragma once
#include "Solver.h"

using namespace utils;

struct KurganovFluxMethod : BaseFluxMethod {
protected:
  using BaseFluxMethod::BaseFluxMethod;
  Array F(index e, index from, index to);
};