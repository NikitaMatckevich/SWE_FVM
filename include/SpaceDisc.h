#pragma once
#include "MUSCLObject.h"
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

	SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, const VolumeField& v0)
    : MUSCLObject(std::move(b), v0)
    , m_fluxer(fluxer)
	{
		m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
	}

	SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, VolumeField&& v0)
    : MUSCLObject(std::move(b), std::move(v0))
    , m_fluxer(fluxer)
  {
		m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
	}

	inline double GetMinLenToWavespeed()  const { return m_min_length_to_wavespeed; }	
	inline const  Storage<3>& GetFluxes() const { return m_f; }
	
	void ComputeFluxes();
	void DumpFluxes(const std::string& filename) const;
	
 protected:
	
	Storage<3> m_f;
  Fluxer 		 m_fluxer;

	double       m_min_length_to_wavespeed;
  const double m_min_wavespeed_to_capture = 1e-10;
};
