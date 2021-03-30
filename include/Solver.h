#pragma once
#include <MUSCLObject.h>
#include <memory>
#include <iostream>

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

Array ElemFlux(const Eigen::Vector2d& n, const Array& Ucons);

struct BaseSpaceDisc : MUSCLObject {

	BaseSpaceDisc(Bathymetry&& b, const VolumeField& v0, double cor = 0., double tau = 0.)
    : MUSCLObject(std::move(b), v0)
		, cor_(cor)
		, tau_(tau)
  {
		F_.resize(Eigen::NoChange, this->mesh().num_edges());
	}

	BaseSpaceDisc(Bathymetry&& b, VolumeField&& v0, double cor = 0., double tau = 0.)
    : MUSCLObject(std::move(b), std::move(v0))
		, cor_(cor)
		, tau_(tau)
  {
		F_.resize(Eigen::NoChange, this->mesh().num_edges());
	}

	inline double getMinLenToWavespeed() const { return min_length_to_wavespeed_; }
	Array RHS(idx i);

 protected:

	virtual Array UserFlux(idx e, idx from, idx to) = 0;
	Array Flux(idx e, idx from, idx to);
  	
	Storage3d F_;
	double       min_length_to_wavespeed_ = 1.;
  const double min_wavespeed_to_capture_ = 1e-10;
	double cor_;                        // Coriolis force parameter
  double tau_;                        // bottom friction parameter
};

struct BaseTimeSolver {

	BaseTimeSolver(std::unique_ptr<BaseSpaceDisc> SpDisc)
		: SpDisc_(std::move(SpDisc)) {};

	double CFLdt() const {
    return CFLconst_ * SpDisc_->getMinLenToWavespeed();
  }
	void Step(double dt) {
		SpDisc_->reconstruct_all();
		UserStep(dt);
	}
	virtual void UserStep(double dt) = 0;
	BaseSpaceDisc* getSpaceDisc() { return SpDisc_.get(); }

 protected:
	std::unique_ptr<BaseSpaceDisc> SpDisc_;
	const double CFLconst_ = 0.9;   // time step is 90% of CFL
};
