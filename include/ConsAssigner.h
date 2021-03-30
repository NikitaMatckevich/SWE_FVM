#pragma once
#include <Bathymetry.h>

struct ConsAssigner {
	using Array = Eigen::Array3d;
	using Field = Eigen::Array3Xd;

	ConsAssigner(Field* field, double b, idx col);
	const Array& operator =(const Array& rhs);
	const Array& operator+=(const Array& rhs);
	const Array& operator-=(const Array& rhs);
	const Array& operator*=(double scalar);
	const Array& operator/=(double scalar);

 private:
	Field* const field_;
	double b_;
	idx col_;
};
