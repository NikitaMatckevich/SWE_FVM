#pragma once
#include "SpaceDisc.h"

struct TimeDisc {
	
	using Solver = std::function<void(TimeDisc * const, double)>;

	TimeDisc(const Solver& solver, double cor = 0., double tau = 0.)
		: m_solver(solver)
		, m_cor(cor)
		, m_tau(tau) {}

	inline double CFLdt() const { return m_constCFL * m_sd->GetMinLenToWavespeed(); }

	inline SpaceDisc* GetSpaceDisc() { return m_sd; }
	inline const SpaceDisc* GetSpaceDisc() const { return m_sd; }
  inline void SetSpaceDisc(SpaceDisc* sd) { m_sd = sd; }

	Array<3> RHS(Idx i, double dt) const;

	inline void Step(double dt) {
    if (m_sd) {
		  m_sd->ReconstructAll();
		  m_sd->ComputeFluxes();
		  m_solver(this, dt);
    }
	}

	inline const double GetTau() const { return m_tau; }
	inline const double GetCor() const { return m_cor; }

 protected:

	double ComputeDrainingDt(Idx i) const;
 
	SpaceDisc*   m_sd {nullptr};
	Solver 		   m_solver;
	const double m_constCFL = 0.15;  // time step coef. is < 1/6 [Liu et al., 2018]
	double 			 m_cor;              // rotation force parameter
  double 			 m_tau;              // bottom friction parameter
};
