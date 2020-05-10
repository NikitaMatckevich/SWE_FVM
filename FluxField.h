#pragma once
#include "MUSCLObject.h"

namespace {
  using Array = BaseValueField::Array;

}

template <class FluxMethod>
struct FluxField {
  using Array = BaseValueField::Array;
  MUSCLObject U;
  Array Flux(Eigen::Vector2d n, Array Ucons) {
    Array res(Array::Zero());
    double h = Ucons[0];
    if (is_wet(h)) {
      auto hvel = Ucons.tail<2>();
      double hveln = vel.matrix().dot(n);
      res[0] = hveln;
      res.tail<2>() = (hveln / h) * hvel  + (0.5 * h * h) * n.array();
    }
    return res;
  }
  Array F() {
    Array res(Array::Zero());
    auto m = U.b.mesh();
    index t0 = m.edge_triangs(e)[0];
    index t1 = m.edge_triangs(e)[1];
    auto n = m.normal(e, );
    if (m.is_edge_boundary(e)) {
      switch (t1) {
      case SOLID_WALL: res = Flux(n, );
      }
    }
    else if (dry_cell(t0) && dry_cell(t1)) {
      return res;
    }
    else {
      Ul = U.Edg.prim(e, t0, t1);
      Ur = U.Edg.prim(e, t1, t0);
      auto is_tvd_ok = abs(Ur - Ul) < abs(U.Vol.prim(t1) - U.Vol.prim(t0));
      Ul = is_tvd_ok.select(Ul, U.Vol.prim(t0));
      Ur = is_tvd_ok.select(Ur, U.Vol.prim(t1));
      double cl = sqrt(Ul[0] - b.e(e));
      double cr = sqrt(Ur[0] - b.e(e));
      double laml_in  = Ul.tail<2>().matrix().dot(n) - cl;
      double lamr_in = Ur.tail<2>().matrix().dot(n) - cr;
      double a_in = 

      double a_out(Vector2d const& n, double B,
        Array3d const& Y1, Array3d const& Y2) {
        double c1 = sqrt(h(Y1[0], B));
        double c2 = sqrt(h(Y2[0], B));
        double lam_in = u(Y1).dot(n) + c1;
        double lam_out = u(Y2).dot(n) + c2;
        return std::max<double>({ lam_in, lam_out, 0. });
      }

      double ak_in = a_in(n, B_.e(e), Y0, Y1);
      double ak_out = a_out(n, B_.e(e), Y0, Y1);
      if (ak_in + ak_out > min_wavespeed_to_capture_) {
        ///////////////////////////////////////////////////
        double r0 = 2. * m.area(t0) / m.l(e);
        double r1 = 2. * m.area(t1) / m.l(e);
        double length_to_wavespeed =
          std::min(r0, r1) / std::max(ak_in, ak_out);
        if (length_to_wavespeed < min_length_to_wavespeed_)
          min_length_to_wavespeed_ = length_to_wavespeed;
        ///////////////////////////////////////////////////
        F_.col(e) = (ak_in * Flux(n, Edg_.cons.) +
          ak_out * Flux(n, B_.e(e), Y0) +
          ak_in * ak_out * (::U(Y0, B_.e(e)) - ::U(Y1, B_.e(e))))
          / (ak_in + ak_out);
      }
    }
    return res;
  }

  double min_length_to_wavespeed() const { return min_length_to_wawespeed_; }

private:
  Bathymetry const& B_;
  VolumeField const& Vol_;
  EdgeField const& Edg_;
  double min_length_to_wawespeed_;
};