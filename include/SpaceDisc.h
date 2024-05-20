#pragma once
#include <MUSCLObject.h>
#include <Tests.h>
#include <fstream>
#include <iostream>
#include <string>

inline Array<3> ElemFlux(const Eigen::Vector2d& n, const Array<3>& U) {

  Array<3> res = Array<3>::Zero();

  double h = U[0];

  if (IsWet(h)) {
    Array<2> hvel = U.tail<2>();
    double hveln = hvel.matrix().dot(n);
    res[0] = hveln;
    res.tail<2>() = (hveln / h) * hvel + (0.5 * h * h) * n.array();
  }

  return res;
}

struct SpaceDisc : MUSCLObject {

    using Fluxer = std::function<Array<3>(SpaceDisc* const, Idx, Idx, Idx, double*)>;

	SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, const VolumeField& v0, double cor = 0, double tau = 0);
	SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, VolumeField&& v0, double cor = 0, double tau = 0);
  
	inline const EdgeField& GetEdgField() const noexcept { return m_edg; }
	inline const EdgeField& GetSrcField() const noexcept { return m_src; }
    inline const Storage<3>& GetFluxes() const noexcept { return m_f; }
    inline const double GetTau() const noexcept { return m_tau; }
	inline const double GetCor() const noexcept { return m_cor; }
    inline const double GetMinLenToWavespeed() const noexcept { return m_min_length_to_wavespeed; }	

    void ComputeInterfaceValues();
	void ComputeFluxes();

 protected:	
    void UpdateInterfaceValues(const MUSCL& muscl);

    EdgeField   m_edg;
	EdgeField   m_src;
    Storage<3>  m_f;
    Fluxer 		  m_fluxer;
  
    double       m_min_length_to_wavespeed;
    const double m_min_wavespeed_to_capture = 1e-10;
    const double m_cor; // rotational force parameter
    const double m_tau; // bottom friction parameter
};
