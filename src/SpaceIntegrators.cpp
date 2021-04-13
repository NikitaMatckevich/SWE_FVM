#include <SpaceIntegrators.h>

Array<3> KurganovSpaceDisc::UserFlux(idx e, idx from, idx to) {

	const auto& m = mesh();
	const auto& Edg = Edg_;	

	auto n = m.norm(e, from);
 
	double ul = Edg.prim(e, from, to).tail<2>().matrix().dot(n);
	double cl = sqrt(Edg.h(e, from, to));

	double ur = Edg.prim(e, to, from).tail<2>().matrix().dot(n);
	double cr = sqrt(Edg.h(e, to, from));

	double al = -std::min({ ul - cl, ur - cr, 0. });
	double ar =  std::max({ ul + cl, ur + cr, 0. });

  // std::cout << from << ' ' << to << ' ' << ' ' << al << ' ' << ar << '\n';
	
	if (al + ar <= min_wavespeed_to_capture_) {
		return Array<3>::Zero();
	} 

	double dl = 2. * m.area(from) / m.l(e);
	double dr = 2. * m.area(to)   / m.l(e);
	double length_to_wavespeed = std::min(dl, dr) / std::max(al, ar);
	if (length_to_wavespeed < min_length_to_wavespeed_) {
		min_length_to_wavespeed_ = length_to_wavespeed;
	}

	/* TODO: delete */
  // al = ar = 0.5;
	/**/

	const auto& Ul = Edg.cons(e, from, to);
	const auto& Ur = Edg.cons(e, to, from);

	//return ElemFlux(n, Ul);

	return (al * ElemFlux(n, Ur) + ar * ElemFlux(n, Ul) - (al * ar) * (Ur - Ul)) / (al + ar);
	//return (al * ElemFlux(n, Ur) + ar * ElemFlux(n, Ul)) / (al + ar);
}
