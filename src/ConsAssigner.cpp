#include <ConsAssigner.h>

Array<3> PrimAssigner::get() const {
	return std::as_const(*field_).col(col_);
}

void PrimAssigner::operator =(const Array<3>& rhs) {
	double h = rhs(0) - b_;
	field_->col(col_) = rhs;
	if (h < eta_1) {
		field_->col(col_).tail<2>() *= (sqrt(2) * h / sqrt(h*h + eta_2));
	}
}

Array<3> ConsAssigner::get() const {
	Array<3> res = std::as_const(*field_).col(col_);
	res.tail(2) *= (res(0) -= b_);
	return res;
}

void ConsAssigner::operator =(const Array<3>& rhs) {
	double h = rhs(0);
	double ih;
	if (h < eta_1) {
		ih = sqrt(2) * h / sqrt(h*h*h*h + eta_4);
	} else {
		ih = 1. / h;
	}
	field_->col(col_) = Array<3>{ h + b_, rhs(1) * ih, rhs(2) * ih };
}
