#pragma once
#include <MUSCLObject.h>
#include <fstream>
#include <iostream>
#include <string>

inline Array<3> ElemFlux(const Eigen::Vector2d& n, const Array<3>& Ucons) {
  Array<3> res = Array<3>::Zero();
  double h = Ucons[0];
  if (is_wet(h)) {
    Array<2> hvel = Ucons.tail<2>();
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
		F_.resize(Eigen::NoChange, this->mesh().num_edges());
	}

	BaseSpaceDisc(Bathymetry&& b, VolumeField&& v0)
    : MUSCLObject(std::move(b), std::move(v0))
  {
		F_.resize(Eigen::NoChange, this->mesh().num_edges());
	}

	double getMinLenToWavespeed() const { return min_length_to_wavespeed_; }	
	const  Storage<3>& getFluxes() const { return F_; }
	
	void compute_fluxes() {

		min_length_to_wavespeed_ = 1.;
		auto const& m = mesh();

		for (size_t i = 0; i < m.num_edges(); ++i) {

			auto const& it = m.edge_triangs(i);
			idx lf = it[0];
			idx lt = it[1];

			auto n = m.norm(i, lf);
			
			switch (lt) {
			case static_cast<idx>(boundaries::SOLID_WALL):
				F_.col(i) = ElemFlux(n, Array<3>{ Vol_.h(lf), 0., 0. });
				break;
			default:
				F_.col(i) = static_cast<Derived&>(*this).UserFlux(i, lf, lt);
			}
		}
	}

	void dump(const std::string& filename) const {
		std::ofstream fout(filename);
		fout << "x\ty\tf1\tf2\tf3\tt1\tt2\n";
		for (size_t i = 0; i < mesh().num_edges(); i++) {
			const auto& e = mesh().e(i);
			const auto& it = mesh().edge_triangs(i);
			fout << e[0] << '\t' << e[1] << '\t';
			fout << F_(0, i) << '\t' << F_(1, i) << '\t' << F_(2, i) << '\t';
			fout << it[0] << '\t' << it[1] << '\n';
		}
		fout.close();
	}
	
 protected:
	
	Storage<3> F_;

	double       min_length_to_wavespeed_;
  const double min_wavespeed_to_capture_ = 1e-10;
};
