#include "ValueFieldAssigner.h"
#include "Solver.h"

using Array = Assigner::Array;
using Field = Assigner::Field;

Assigner::Assigner(Field* field, const Bathymetry& b, index col)
	: field_(field), b_(b), col_(col) {}
const Array& Assigner::operator =(const Array& rhs) {
	double h = rhs(0);
	double b = b_.t(col_);
	field_->col(col_) = is_wet(h) ? Array{ h + b, rhs(1) / h, rhs(2) / h } : dry_state(b);
	return rhs;
}
const Array& Assigner::operator+=(const Array& rhs) {
	return this->operator=(std::as_const(*field_).col(col_) + rhs);
}
const Array& Assigner::operator-=(const Array& rhs) {
	return this->operator=(std::as_const(*field_).col(col_) - rhs);
}
const Array& Assigner::operator*=(double scalar) {
	return this->operator=(std::as_const(*field_).col(col_) * scalar);
}
const Array& Assigner::operator/=(double scalar) {
	if (scalar == 0.)
		throw SolverError("Division by zero when writing in value field in conservative form");
	return this->operator*=(1. / scalar);
}