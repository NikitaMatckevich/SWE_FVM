#include "KurganovFluxMethod.h"
#include <algorithm>

using Array = BaseFluxMethod::Array;

Array KurganovFluxMethod::F(index e, index from, index to) {
  Array res(Array::Zero());
  auto const& m = bathymetry().mesh();
  auto n = m.norm(e, from);
  if (m.is_edge_boundary(e)) {
    switch (to) {
    case SOLID_WALL: res = ElemFlux(n, Array{ Vol_.h(from), 0., 0. });
    }
  }
  else if (is_dry_cell(from) && is_dry_cell(to)) {
    return res;
  }
  else {
    auto Vtmin = Vol_.prim(to).min(Vol_.prim(from));
    auto Vtmax = Vol_.prim(to).max(Vol_.prim(from));
    auto TVD = [&Vtmin, &Vtmax](const auto& Vnew, const auto& Vold)->Array {
      return ((Vtmin <= Vnew) && (Vnew <= Vtmax)).select(Vnew, Vold);
    };
    Array Vl = TVD(Edg_.prim(e, from, to), Vol_.prim(from));
    Array Vr = TVD(Edg_.prim(e, to, from), Vol_.prim(to));
    double ul = Vl.tail<2>().matrix().dot(n);
    double cl = sqrt(Vl[0] - bathymetry().e(e));
    double ur = Vr.tail<2>().matrix().dot(n);
    double cr = sqrt(Vr[0] - bathymetry().e(e));
    double al = -std::min({ ul - cl, ur - cr, 0. });
    double ar = std::max({ ul + cl, ur + cr, 0. });

    Edg_.prim(e, from, to) = Vl;
    Edg_.prim(e, to, from) = Vr;

    if (al + ar > min_wavespeed_to_capture_) {
      double dl = 2. * m.area(from) / m.l(e);
      double dr = 2. * m.area(to) / m.l(e);
      double length_to_wavespeed = std::min(dl, dr) / std::max(al, ar);
      if (length_to_wavespeed < min_length_to_wavespeed_)
        min_length_to_wavespeed_ = length_to_wavespeed;
      const auto& Ul = std::as_const(Edg_).cons(e, from, to);
      const auto& Ur = std::as_const(Edg_).cons(e, to, from);
      res = al * ElemFlux(n, Ur) + ar * ElemFlux(n, Ul) + al * ar / (al + ar) * (Ur - Ul);
    }
  }
  return res;
}