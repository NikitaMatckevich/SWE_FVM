#pragma once
#include "SpaceDisc.h"

struct TimeDisc {
	
	explicit TimeDisc(SpaceDisc * sd = nullptr)
    : m_sd(sd) {}	
  
  inline       SpaceDisc* GetSpaceDisc()       { return m_sd; }
	inline const SpaceDisc* GetSpaceDisc() const { return m_sd; }
  inline void             SetSpaceDisc(SpaceDisc* sd) { m_sd = sd; }
	
  inline double CFLdt() const { return m_constCFL * m_sd->GetMinLenToWavespeed(); }
	
	Array<3> RHS(Idx i, double dt) const;

 protected:

	double ComputeDrainingDt(Idx i) const;
 
	SpaceDisc*   m_sd;
	const double m_constCFL = 0.15;  // time step coef. is < 1/6 [Liu et al., 2018]
};
