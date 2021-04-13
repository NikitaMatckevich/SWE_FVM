#pragma once
#include <SpaceDisc.h>

template <class SpaceDisc>
struct TimeDisc {
	
	static_assert(std::is_convertible_v<SpaceDisc*, BaseSpaceDisc<SpaceDisc>*>);
	using solver_t = std::function<void(TimeDisc<SpaceDisc> * const, double)>;

	TimeDisc(SpaceDisc* sd, const solver_t& solver, double cor = 0., double tau = 0.)
		: sd(sd)
		, solver(solver)
		, cor_(cor)
		, tau_(tau) {}

	double CFLdt() const { return CFLconst_ * sd->getMinLenToWavespeed(); }

	SpaceDisc* getSpaceDisc() { return sd; }
	const SpaceDisc* getSpaceDisc() const { return sd; }

	Array<3> RHS(idx i, double dt) const {

		Array<3> res = { 0., cor_ * sd->getVolField().v(i), -cor_ * sd->getVolField().u(i) };

		auto const& m = sd->mesh();
		auto const& ie = m.triang_edges(i);
		auto const& it = m.triang_triangs(i);
		auto const& F = sd->getFluxes();

		double iS = 1. / m.area(i);

		double dti = compute_draining_dt(i);
		
		for (int k = 0; k < 3; k++) {

			int sgn = m.edge_triangs(ie[k])[0] == i ? 1 : -1;
			
			double dtik = compute_draining_dt(it[k]);
			double dtk = (sgn * F(0, ie[k])) > 0. ? std::min(dt, dti) : std::min(dt, dtik);
	
			double c_ek = iS * m.l(ie[k]);
 			double h_ek = sd->getEdgField().h(ie[k], i, it[k]);
			
			res -= dtk * sgn * c_ek * F.col(ie[k]);
			//res.tail<2>() += dtk * (m.norm(ie[k], i) * c_ek * (0.5 * h_ek * h_ek)).array();
		}
		
		res.tail<2>() -= dt * sd->bath().grad(i) * sd->getVolField().h(i);

		return res;
	}

	void Step(double dt) {
		sd->reconstruct_all();
		sd->compute_fluxes();
		solver(this, dt);
	}

 protected:

	double compute_draining_dt(idx i) const {

		if (i < 0)
			return std::numeric_limits<double>::infinity(); 

		if (sd->is_dry_cell(i))
			return 0.;

		auto const& m = sd->mesh();
		auto const& ie = m.triang_edges(i);
		auto const& F = sd->getFluxes();

    double sum = 0.;
    for (int k = 0; k < 3; k++) {
			idx itk = m.edge_triangs(ie[k])[0];
			double f_ek = F(0, ie[k]);
      sum += std::max(0., (i == itk ? f_ek : -f_ek));
    }

		return sum > 0. ? 
			m.area(i) * sd->getVolField().h(i) / sum :
			std::numeric_limits<double>::infinity()  ;
	}

	SpaceDisc* sd;
	solver_t solver;
	const double CFLconst_ = 0.15;   // time step is 90% of CFL

	double cor_;                     // Coriolis force parameter
  double tau_;                     // bottom friction parameter
};
