#pragma once
#include "SpaceDisc.h"

template <class SpaceDisc>
struct TimeDisc {
	
	static_assert(std::is_convertible_v<SpaceDisc*, BaseSpaceDisc<SpaceDisc>*>);
	using Solver = std::function<void(TimeDisc<SpaceDisc>* const, double)>;

	TimeDisc(SpaceDisc* sd, const Solver& solver, double cor = 0., double tau = 0.)
		: m_sd(sd)
		, m_solver(solver)
		, m_cor(cor)
		, m_tau(tau) {}

	double CFLdt() const { return m_constCFL * m_sd->GetMinLenToWavespeed(); }

	SpaceDisc* GetSpaceDisc() { return m_sd; }
	const SpaceDisc* GetSpaceDisc() const { return m_sd; }

	Array<3> RHS(Idx i, double dt) const {

		const auto& vol = m_sd->GetVolField();
		const auto& s   = m_sd->GetSrcField();
		const auto& f   = m_sd->GetFluxes();

		Array<3> res = { 0., m_cor * vol.v(i), -m_cor * vol.u(i) };

		const auto& m  = m_sd->Mesh();
		const auto& ie = m.TriangEdges(i);
		const auto& it = m.TriangTriangs(i);
		
		double i_area = 1. / m.Area(i);

		double dti = ComputeDrainingDt(i);
		Array<3> dts = (dt / 3.) * s.prim(i); 
	
		for (int k = 0; k < 3; k++) {

			int sgn = m.EdgeTriangs(ie[k])[0] == i ? 1 : -1;
			
			double dtik = ComputeDrainingDt(it[k]);
			double dtk = (sgn * f(0, ie[k])) > 0. ? std::min(dt, dti) : std::min(dt, dtik);
	
			double c_ek = i_area * m.L(ie[k]);
 			double h_ek = m_sd->GetEdgField().h(ie[k], i, it[k]);
			
			res -= dtk * sgn * c_ek * f.col(ie[k]);

			// TODO(nikitamatckevich): if doesn't work, delete lines 51, uncomment line 54
			//res.tail<2>() += dtk * (m.Norm(ie[k], i) * c_ek * (0.5 * h_ek * h_ek)).array();
		}
		
		//res -= dt * s.prim(i) * vol.h(i);
		res.tail<2>() -= dt * m_sd->Bath().Gradient(i) * vol.h(i);

		return res;
	}

	void Step(double dt) {
		m_sd->ReconstructAll();
		m_sd->ComputeFluxes();
		m_solver(this, dt);
	}

	const double GetTau() const { return m_tau; }
	const double GetCor() const { return m_cor; }

 protected:

	double ComputeDrainingDt(Idx i) const {

		if (i < 0)
			return std::numeric_limits<double>::infinity(); 

		if (m_sd->IsDryCell(i))
			return 0.;

		const auto& m  = m_sd->Mesh();
		const auto& ie = m.TriangEdges(i);
	
		const auto& f = m_sd->GetFluxes();

    double sum = 0.;
    for (int k = 0; k < 3; k++) {
			Idx itk = m.EdgeTriangs(ie[k])[0];
			double f_ek = f(0, ie[k]);
      sum += std::max(0., (i == itk ? f_ek : -f_ek));
    }

		return sum > 0. ? 
			m.Area(i) * m_sd->GetVolField().h(i) / sum :
			std::numeric_limits<double>::infinity();
	}

	SpaceDisc*   m_sd;
	Solver 		   m_solver;
	const double m_constCFL = 0.15;  // time step coef. is < 1/6 [Liu et al., 2018]
	double 			 m_cor;              // rotation force parameter
  double 			 m_tau;              // bottom friction parameter
};
