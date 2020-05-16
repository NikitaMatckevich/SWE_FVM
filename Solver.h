#pragma once
#include "MUSCLObject.h"

MUSCLObject::Array ElemFlux(const Eigen::Vector2d& n, const MUSCLObject::Array& Ucons);

template <class DerivedFluxMethod>
struct BaseFluxMethod : MUSCLObject {
  using Derived = DerivedFluxMethod;
protected:
  using MUSCLObject::MUSCLObject;

  Array F(index e, index from, index to) {
    Array res(Array::Zero());
    if (mesh().is_edge_boundary(e)) {
      auto const& m = mesh();
      auto n = m.norm(e, from);
      switch (to) {
      case SOLID_WALL: res = ElemFlux(n, Array{ Vol_.h(from), 0., 0. });
      }
    }
    else if (!is_dry_cell(from) || !is_dry_cell(to))
      res = static_cast<Derived*>(this)->InnerFlux(e, from, to);
    return res;
  }

  double       min_length_to_wavespeed_ = 1.;
  const double min_wavespeed_to_capture_ = 1e-10;
};

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

template <class FluxMethod>
struct BaseSolver : FluxMethod {
  using Array = BaseFluxMethod::Array;
  BaseSolver(Bathymetry&& b, VolumeField&& v0, double cor = 0., double tau = 0.)
    : FluxMethod(std::move(b), std::move(v0))
    , cor_(cor)
    , tau_(tau)
  {}
  double CFL_dt() const {
    return CFL_constant_ * min_length_to_wavespeed_;
  }
  Array rhs(index i) {
    Array res = { 0., cor_ * Vol_.v(i), -cor_ * Vol_.u(i) };
    auto const& m = mesh();
    auto const& ie = m.triang_edges(i);
    auto const& it = m.triang_triangs(i);
    double iS = 1. / m.area(i);
    for (int k = 0; k < 3; k++)
      res -= iS * m.l(ie[k]) * F(ie[k], i, it[k]);
    return res;
  }
private:
  const double CFL_constant_ = 0.9;   // time step is 90% of CFL
  double cor_;                        // Coriolis force parameter
  double tau_;                        // bottom friction parameter
};