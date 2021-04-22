#pragma once
#include <SpaceDisc.h>

namespace Wavespeeds {
  Array<2> Rusanov(double ul, double hl, double ur, double hr); 
  Array<2> Davis(double ul, double hl, double ur, double hr); 
  Array<2> Einfeldt(double ul, double hl, double ur, double hr);
}

template <Array<2> (*WavespeedCalc)(double, double, double, double)>
Array<3> HLLFlux(SpaceDisc * const sd, Idx e, Idx from, Idx to, double * r) {

	const auto& m = sd->Mesh();
	const auto& edg = sd->GetEdgField();	

	auto n = m.Norm(e, from);
 
	double ul = edg.vel(e, from, to).matrix().dot(n);
	double hl = edg.h(e, from, to);

	double ur = edg.vel(e, to, from).matrix().dot(n);
	double hr = edg.h(e, to, from);

  constexpr double min_depth_to_capture = 1e-10;
	if (hl + hr <= min_depth_to_capture) {
		return Array<3>::Zero();
	}

  Array<2> a = WavespeedCalc(ul, hl, ur, hr);
	double al = std::min(0., a[0]);
	double ar = std::max(0., a[1]);

  constexpr double min_wavespeed_to_capture = 1e-10;
	if (ar - al <= min_wavespeed_to_capture) {
		return Array<3>::Zero();
	} 

  if (r) {
    double dl = 2. * m.Area(from) / m.L(e);
	  double dr = 2. * m.Area(to)   / m.L(e);
	  double length_to_wavespeed = std::min(dl, dr) / std::max(-al, ar);
		*r = std::min(*r, length_to_wavespeed);
  }

	const auto& Ul = edg.cons(e, from, to);
	const auto& Ur = edg.cons(e, to, from);

	return (ar * ElemFlux(n, Ul) - al * ElemFlux(n, Ur) + (al * ar) * (Ur - Ul)) / (ar - al);
}

template <Array<2> (*WavespeedCalc)(double, double, double, double)>
Array<3> HLLCFlux(SpaceDisc * const sd, Idx e, Idx from, Idx to, double * r) {

  const auto& m = sd->Mesh();
	const auto& edg = sd->GetEdgField();	

  auto t = m.Tang(e, from);
	auto n = m.Norm(e, from);
  
	double ul = edg.vel(e, from, to).matrix().dot(n);
  double vl = edg.vel(e, from, to).matrix().dot(t);
  double hl = edg.h(e, from, to);

	double ur = edg.vel(e, to, from).matrix().dot(n);
	double vr = edg.vel(e, to, from).matrix().dot(t);
  double hr = edg.h(e, to, from);

  constexpr double min_depth_to_capture = 1e-10;
	if (hl + hr <= min_depth_to_capture) {
		return Array<3>::Zero();
	}

  Array<2> a = WavespeedCalc(ul, hl, ur, hr);
	double al = a[0];
	double ar = a[1];

  double ustar = (ar - ur) * hr * ur - (al - ul) * hl * ul + 0.5 * (hl * hl - hr * hr);
  ustar /= (hr * (ar - ur) - hl * (al - ul));

  if (r) {
    double dl = 2. * m.Area(from) / m.L(e);
	  double dr = 2. * m.Area(to)   / m.L(e);
	  double length_to_wavespeed = std::min(dl, dr) / std::max(tol, std::max(al, ar));
		*r = std::min(*r, length_to_wavespeed);
  }

	const auto& Ul = edg.cons(e, from, to);
	const auto& Ur = edg.cons(e, to, from);

  if (ustar <= 0) {
    Array<2> vel = Array<2>{vr, ustar};
    double urstar = vel.matrix().dot(t);
    double vrstar = vel.matrix().dot(n);
    double hrstar = hr * (ar - ur) / (ar - ustar);
    Array<3> Urstar = Array<3>{hrstar, hrstar * urstar, hrstar * vrstar};
    return ElemFlux(n, Ur) + std::max(0., ar) * (Urstar - Ur);
  }
  else {
    Array<2> vel = Array<2>{vl, ustar};
    double ulstar = vel.matrix().dot(t);
    double vlstar = vel.matrix().dot(n);
    double hlstar = hl * (al - ul) / (al - ustar);
    Array<3> Ulstar = Array<3>{hlstar, hlstar * ulstar, hlstar * vlstar};
    return ElemFlux(n, Ul) + std::min(0., al) * (Ulstar - Ul);
  }
}
