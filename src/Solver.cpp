#include <Solver.h>

Array ElemFlux(const Eigen::Vector2d& n, const Array& Ucons) {
  Array res(Array::Zero());
  double h = Ucons[0];
  if (is_wet(h)) {
    auto hvel = Ucons.tail<2>();
    double hveln = hvel.matrix().dot(n);
    res[0] = hveln;
    res.tail<2>() = (hveln / h) * hvel + (0.5 * h * h) * n.array();
  }
  return res;
}

Array BaseSpaceDisc::Flux(idx e, idx from, idx to) {
	auto const& m = mesh();
	auto n = m.norm(e, from);
	if (m.is_edge_boundary(e)) {
		switch (to) {
		case static_cast<idx>(boundaries::SOLID_WALL):
			F_.col(e) = ElemFlux(n, Array{ Vol_.h(from), 0., 0. });
			break;
		}
	}
	else if (is_dry_cell(from) && is_dry_cell(to)) {
		F_.col(e) = Array::Zero();
	}
	else {
		F_.col(e) = UserFlux(e, from ,to);
	}
	return F_.col(e);
}

Array BaseSpaceDisc::RHS(idx i) {
	Array res = { 0., cor_ * Vol_.v(i), -cor_ * Vol_.u(i) };
	auto const& m = mesh();
	auto const& ie = m.triang_edges(i);
	auto const& it = m.triang_triangs(i);
	double iS = 1. / m.area(i);
	for (int k = 0; k < 3; k++)
		res -= iS * m.l(ie[k]) * Flux(ie[k], i, it[k]);
	return res;
}
