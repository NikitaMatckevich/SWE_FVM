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

template <class Derived>
struct BaseSpaceDisc : MUSCLObject {

	BaseSpaceDisc(Bathymetry&& b, const VolumeField& v0)
    : MUSCLObject(std::move(b), v0)
	{
		m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
	}

	BaseSpaceDisc(Bathymetry&& b, VolumeField&& v0)
    : MUSCLObject(std::move(b), std::move(v0))
  {
		m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
	}

	double GetMinLenToWavespeed()  const { return m_min_length_to_wavespeed; }	
	const  Storage<3>& GetFluxes() const { return m_f; }
	
	void ComputeFluxes() {

		m_min_length_to_wavespeed = 1.;
		const auto& m = Mesh();

		for (size_t i = 0; i < m.NumEdges(); ++i) {

			const auto& it = m.EdgeTriangs(i);
			Idx lf = it[0];
			Idx lt = it[1];

			auto n = m.Norm(i, lf);
			
			switch (lt) {
			case static_cast<Idx>(Boundaries::SOLID_WALL):
				m_f.col(i) = ElemFlux(n, Array<3>{ m_vol.h(lf), 0., 0. });
				break;
			default:
				m_f.col(i) = static_cast<Derived&>(*this).UserFlux(i, lf, lt);
			}
		}
	}

	void DumpFluxes(const std::string& filename) const {

		std::ofstream fout(filename);

		fout << "x\ty\tf1\tf2\tf3\tt1\tt2\n";

		for (size_t i = 0; i < Mesh().NumEdges(); i++) {
			
      const auto& e = Mesh().E(i);
			const auto& it = Mesh().EdgeTriangs(i);
			
      fout << e[0] << '\t' << e[1] << '\t';
			fout << m_f(0, i) << '\t' << m_f(1, i) << '\t' << m_f(2, i) << '\t';
			fout << it[0] << '\t' << it[1] << '\n';
		}

		fout.close();
	}
	
 protected:
	
	Storage<3> m_f;

	double       m_min_length_to_wavespeed;
  const double m_min_wavespeed_to_capture = 1e-10;
};
