#pragma once
#include "Bathymetry.h"

template <class Derived, size_t k>
struct BaseAssigner {
	
	BaseAssigner(Storage<k>* field, double b, Idx col)
		: m_field(field), m_b(b), m_col(col) {}
	
	void operator+=(const Array<k>& rhs) {
		Child().operator=(Child().Get()  + rhs);
	}
	void operator-=(const Array<k>& rhs) {
		Child().operator=(Child().Get() - rhs);
	}
	void operator*=(double scalar) {
		Child().operator=(Child().Get() * scalar);
	}
	void operator/=(double scalar) {
		if (scalar == 0.) {
			throw std::runtime_error("Division by zero when writing in value field in conservative form");
    }
		Child().operator*=(1. / scalar);
	}

 protected:
	Storage<k>* const m_field;
	double m_b;
	Idx m_col;

	const double m_eta_1 = 1e-3;
	const double m_eta_2 = 1e-6;
	const double m_eta_4 = 1e-12;

 private:	
  inline Derived& Child() { return static_cast<Derived&>(*this); }
};

struct PrimAssigner : BaseAssigner<PrimAssigner, 3> {
	using BaseAssigner<PrimAssigner, 3>::BaseAssigner;
	Array<3> Get() const;
	void operator=(const Array<3>& rhs);
};

struct ConsAssigner : BaseAssigner<ConsAssigner, 3> {
	using BaseAssigner<ConsAssigner, 3>::BaseAssigner;
	Array<3> Get() const;
	void operator=(const Array<3>& rhs);
};
