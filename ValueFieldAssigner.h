#pragma once
#include "Bathymetry.h"

struct Assigner {
	using Array = Eigen::Array3d;
	using Field = Eigen::Array3Xd;
	Assigner(Field* field, const Bathymetry& b, index col);
	const Array& operator =(const Array& rhs);
	const Array& operator+=(const Array& rhs);
	const Array& operator-=(const Array& rhs);
	const Array& operator*=(double scalar);
	const Array& operator/=(double scalar);
private:
	Field * const field_;
	Bathymetry const& b_;
	index col_;
};