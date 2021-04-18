#include <Assigners.h>

Array<3> PrimAssigner::Get() const {
	return std::as_const(*m_field).col(m_col);
}

void PrimAssigner::operator =(const Array<3>& rhs) {
	
	double h = rhs(0) - m_b;
	if (!IsWet(h)) {
		m_field->col(m_col) = DryState(m_b);
		return;
	}

	m_field->col(m_col) = rhs;
	if (h < m_eta_1) {
		m_field->col(m_col).tail<2>() *= (sqrt(2) * h / sqrt(h*h + m_eta_2));
	}
}

Array<3> ConsAssigner::Get() const {
	Array<3> res = std::as_const(*m_field).col(m_col);
	res.tail(2) *= (res(0) -= m_b);
	return res;
}

void ConsAssigner::operator =(const Array<3>& rhs) {

	double h = rhs(0);
	if (!IsWet(h)) {
		m_field->col(m_col) = DryState(m_b);
		return;
	}
	double ih;

	if (h < m_eta_1) {
		ih = sqrt(2) * h / sqrt(h*h*h*h + m_eta_4);
	} else {
		ih = 1. / h;
	}

	m_field->col(m_col) = Array<3>{ h + m_b, rhs(1) * ih, rhs(2) * ih };
}
