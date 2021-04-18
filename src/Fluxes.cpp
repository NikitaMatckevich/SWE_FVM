#include <Fluxes.h>

Array<3> KurganovFlux(SpaceDisc * const sd, Idx e, Idx from, Idx to, double * r) {

	const auto& m = sd->Mesh();
	const auto& edg = sd->GetEdgField();	

	auto n = m.Norm(e, from);
 
	double ul = edg.prim(e, from, to).tail<2>().matrix().dot(n);
	double cl = sqrt(edg.h(e, from, to));

	double ur = edg.prim(e, to, from).tail<2>().matrix().dot(n);
	double cr = sqrt(edg.h(e, to, from));

	double al = -std::min({ ul - cl, ur - cr, 0. });
	double ar =  std::max({ ul + cl, ur + cr, 0. });

  // std::cout << from << ' ' << to << ' ' << ' ' << al << ' ' << ar << '\n';
	
  constexpr double min_wavespeed_to_capture = 1e-10;
	if (al + ar <= min_wavespeed_to_capture) {
		return Array<3>::Zero();
	} 

  if (r) {
    double dl = 2. * m.Area(from) / m.L(e);
	  double dr = 2. * m.Area(to)   / m.L(e);
	  double length_to_wavespeed = std::min(dl, dr) / std::max(al, ar);
		*r = std::min(*r, length_to_wavespeed);
  }

	/* TODO: delete */
  // al = ar = 0.5;
	/**/

	const auto& Ul = edg.cons(e, from, to);
	const auto& Ur = edg.cons(e, to, from);

	//return ElemFlux(n, Ul);

	return (al * ElemFlux(n, Ur) + ar * ElemFlux(n, Ul) - (al * ar) * (Ur - Ul)) / (al + ar);
	//return (al * ElemFlux(n, Ur) + ar * ElemFlux(n, Ul)) / (al + ar);
}
