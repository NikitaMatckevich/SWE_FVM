#pragma once
#include "MUSCLObject.h"

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

namespace utils {

  struct BaseFluxMethod : MUSCLObject {
    using Array = BaseValueField::Array;
  protected:
    explicit BaseFluxMethod(VolumeField&& v0);
    double       min_length_to_wavespeed_ = 1.;
    const double min_wavespeed_to_capture_ = 1e-10;
  };

  BaseFluxMethod::Array ElemFlux(const Eigen::Vector2d& n, const BaseFluxMethod::Array& Ucons);

  template <class FluxMethod>
  struct BaseSolver : FluxMethod {
    using Array = BaseFluxMethod::Array;
    BaseSolver(VolumeField&& v0, double cor, double tau)
      : FluxMethod(std::move(v0))
      , cor_(cor)
      , tau_(tau)
    {}
    double CFL_dt() const {
      return CFL_constant_ * min_length_to_wavespeed_;
    }
    Array rhs(index i) {
      Array res = { 0., cor_ * Vol_.v(i), -cor_ * Vol_.u(i) };
      auto const& m = bathymetry().mesh();
      auto const& ie = m.triang_edges(i);
      auto const& it = m.triang_triangs(i);
      double iS = 1. / m.area(i);
      for (int k = 0; k < 3; k++)
        res -= iS * m.l(ie[k]) * F(ie[k], i, it[k]);
      return res;
    }
  private:
    const double CFL_constant_ = 0.9;   // Time step is less than CFL
    double cor_;                        // Coriolis force parameter
    double tau_;                        // Bottom friction parameter
  };
} // namespace utils
