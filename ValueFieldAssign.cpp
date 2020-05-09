#include "ValueField.h"

using Array3d = ValueField::Array3d;

constexpr bool is_wet(double h) noexcept {
	constexpr double h_min = 1e-13; // minimal posible water depth of "wet" cell  
	return h > h_min;
}
Array3d dry_state(double b) noexcept { return { b, 0., 0. }; }

ValueField::AssignColumnCons::AssignColumnCons(ValueField& field, index col) : field_(field), col_(col) {}
const Array3d& ValueField::AssignColumnCons::operator =(const Array3d& rhs) {
	double h = rhs(0);
	double b = field_.b_.t(col_);
	field_.prim(col_) = is_wet(h) ? Array3d{ h + b, rhs(1)/h, rhs(2)/h } : dry_state(b);
	return rhs;
}
const Array3d& ValueField::AssignColumnCons::operator+=(const Array3d& rhs) {
	return this->operator=(std::as_const(field_).cons(col_) += rhs);
}
const Array3d& ValueField::AssignColumnCons::operator-=(const Array3d& rhs) {
	return this->operator=(std::as_const(field_).cons(col_) -= rhs);
}
const Array3d& ValueField::AssignColumnCons::operator*=(double scalar) {
	return this->operator=(std::as_const(field_).cons(col_) *= scalar);
}
const Array3d& ValueField::AssignColumnCons::operator/=(double scalar) {
	if (scalar == 0.)
		throw SolverError("Division by zero when writing in value field in conservative form");
	return this->operator*=(1. / scalar);
}

ValueField::AssignHeightCons::AssignHeightCons(double& target, double b) : target_(target), b_(b) {}
double ValueField::AssignHeightCons::operator =(double h) {
	target_ = h + b_;
	return h;
}
double ValueField::AssignHeightCons::operator+=(double h) {
	return target_ += h;
}
double ValueField::AssignHeightCons::operator-=(double h) {
	return target_ -= h;
}
double ValueField::AssignHeightCons::operator*=(double scalar) {
	return target_ *= scalar;
}
double ValueField::AssignHeightCons::operator/=(double scalar) {
	if (scalar == 0.)
		throw SolverError("Division by zero when writing in value field in conservative form");
	return target_ /= scalar;
}

Array3d ValueField::prim(index i) const { return U_.col(i); }
double  ValueField::w(index i)    const { return U_(0, i); }
double  ValueField::u(index i)    const { return U_(1, i); }
double  ValueField::v(index i)    const { return U_(2, i); }
Array3d ValueField::cons(index i) const { return { h(i), hu(i), hv(i) }; }
double  ValueField::h (index i)   const { return U_(0, i) + b_.t(i); }
double  ValueField::hu(index i)   const { return h(i) * u(i); }
double  ValueField::hv(index i)   const { return h(i) * v(i); }

ValueField::Column3d ValueField::prim(index i) { return U_.col(i); }
double& ValueField::w(index i) { return U_(0, i); }
double& ValueField::u(index i) { return U_(1, i); }
double& ValueField::v(index i) { return U_(2, i); }
ValueField::AssignColumnCons ValueField::cons(index i) { return AssignColumnCons(*this, i); }
ValueField::AssignHeightCons ValueField::h(index i) { return AssignHeightCons(U_(0, i), b_.t(i)); }