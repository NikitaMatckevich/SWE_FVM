#pragma once
#include "SpaceDisc.h"

struct TimeDisc {
	
	explicit TimeDisc(double cor = 0., double tau = 0.)
		: m_cor(cor)
		, m_tau(tau) {}
	
  inline const double GetTau() const { return m_tau; }
	inline const double GetCor() const { return m_cor; }

	inline       SpaceDisc* GetSpaceDisc()       { return m_sd; }
	inline const SpaceDisc* GetSpaceDisc() const { return m_sd; }

  inline void             SetSpaceDisc(SpaceDisc* sd) { m_sd = sd; }
	
  inline double CFLdt() const { return m_constCFL * m_sd->GetMinLenToWavespeed(); }
	
	Array<3> RHS(Idx i, double dt) const;

 protected:

	double ComputeDrainingDt(Idx i) const;
 
	SpaceDisc*   m_sd {nullptr};

	const double m_constCFL = 0.15;  // time step coef. is < 1/6 [Liu et al., 2018]
	double 			 m_cor;              // rotation force parameter
  double 			 m_tau;              // bottom friction parameter
};
