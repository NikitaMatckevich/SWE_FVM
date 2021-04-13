#pragma once
#include <Bathymetry.h>

template <class Derived, size_t k>
struct BaseAssigner {
	
	BaseAssigner(Storage<k>* field, double b, idx col)
		: field_(field), b_(b), col_(col) {}
	
	void operator+=(const Array<k>& rhs) {
		derived().operator=(derived().get()  + rhs);
	}
	void operator-=(const Array<k>& rhs) {
		derived().operator=(derived().get() - rhs);
	}
	void operator*=(double scalar) {
		derived().operator=(derived().get() * scalar);
	}
	void operator/=(double scalar) {
		if (scalar == 0.)
			throw std::runtime_error("Division by zero when writing in value field in conservative form");
		derived().operator*=(1. / scalar);
	}

 protected:
	Storage<k>* const field_;
	double b_;
	idx col_;

	const double eta_1 = 1e-3;
	const double eta_2 = 1e-6;
	const double eta_4 = 1e-12;

 private:
	Derived& derived() { return static_cast<Derived&>(*this); }
};

struct PrimAssigner : BaseAssigner<PrimAssigner, 3> {
	using BaseAssigner<PrimAssigner, 3>::BaseAssigner;
	Array<3> get() const;
	void operator =(const Array<3>& rhs);
};

struct ConsAssigner : BaseAssigner<ConsAssigner, 3> {
	using BaseAssigner<ConsAssigner, 3>::BaseAssigner;
	Array<3> get() const;
	void operator=(const Array<3>& rhs);
};
