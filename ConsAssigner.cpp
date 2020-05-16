#include "ConsAssigner.h"
#include "Solver.h"

using Array = ConsAssigner::Array;
using Field = ConsAssigner::Field;

ConsAssigner::ConsAssigner(Field* field, double b, index col)
	: field_(field), b_(b), col_(col) {}
const Array& ConsAssigner::operator =(const Array& rhs) {
	double h = rhs(0);
	field_->col(col_) = is_wet(h) ? Array{ h + b_, rhs(1) / h, rhs(2) / h } : dry_state(b_);
	return rhs;
}
const Array& ConsAssigner::operator+=(const Array& rhs) {
	return this->operator=(std::as_const(*field_).col(col_) + rhs);
}
const Array& ConsAssigner::operator-=(const Array& rhs) {
	return this->operator=(std::as_const(*field_).col(col_) - rhs);
}
const Array& ConsAssigner::operator*=(double scalar) {
	return this->operator=(std::as_const(*field_).col(col_) * scalar);
}
const Array& ConsAssigner::operator/=(double scalar) {
	if (scalar == 0.)
		throw SolverError("Division by zero when writing in value field in conservative form");
	return this->operator*=(1. / scalar);
}