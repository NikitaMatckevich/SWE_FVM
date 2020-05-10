#pragma once
#include "FluxField.h"

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

template <class ODEMethod, class FluxMethod>
struct Solver {
  double CFL_dt() const {
    return ODEMethod::time_split * CFL_constant_ * F.min_length_to_wavespeed();
  }
  void step(double dt) {
    ODEMethod::step(rhs_, dt);
  }
private:
  VolumeField Uv;
  EdgeField Ue;
  FluxField<FluxMethod> F;
  //RHSObject rhs_;
  const double CFL_constant_ = 0.9;
  double cor_;
  double tau_;
};