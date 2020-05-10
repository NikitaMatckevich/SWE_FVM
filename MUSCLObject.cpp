#include "MUSCLObject.h"
#include "PointOperations.h"
#include "Solver.h"
using namespace Eigen;

namespace {
	Matrix32d gradient_coefs(const Matrix32d& P) {
		Matrix3d L; L << 0,  1, -1,
										-1,  0,  1,
										 1, -1,  0;
		Matrix2d R; R << 0, -1,
										 1,  0;
		Matrix32d res = L * P * R;
		double d = res.topRows<2>().determinant();
		if (abs(d) < tol)
			throw SolverError("Division by zero in gradient coefficient computation");
		res /= d;
		return res;
	}
}

bool MUSCLObject::dry_cell(index i) const {
	return !is_wet(Vol.h(i));
}
bool MUSCLObject::full_wet_cell(index i) const {
	auto const& mesh = b.mesh();
	double max_bp = mesh.triang_points(i).unaryExpr([this](index i) {
		return b.p(i);
	}).maxCoeff();
	return max_bp < Vol.w(i) ? false : mesh.triang_triangs(i).unaryExpr([this](index i) {
		return i >= 0 ? dry_cell(i) : false;
	}).all();
}
bool MUSCLObject::part_wet_cell(index i) const {
	return !dry_cell(i) && !full_wet_cell(i);
}
Matrix32d MUSCLObject::gradient(index i) const {
	assert(full_wet_cell(i));
	auto m = b.mesh();
	auto it = m.triang_triangs(i);
	auto ie = m.triang_edges(i);
	Matrix32d points; // points needed to build reconstruction plane 
	Matrix3d  values; // values at these points
	for (int k = 0; k < 3; k++) {
		if (full_wet_cell(it[k])) {
			points.row(k) = m.t(it[k]);
			values.col(k) = Vol.prim(it[k]).matrix();
		}
		else {
			points.row(k) = m.e(ie[k]);
			values.col(k) = 0.5 * (Vol.prim(i) + Edg.prim(ie[k], it[k], i)).matrix();
		}
	}
	Matrix32d gradient = values * gradient_coefs(points);
	auto const& ip = m.triang_points(i);
	auto dx = (m.p(ip) - m.t(i).replicate<1, 3>()).matrix();
	Array3d wp = Vol.w(i) + (gradient.row(0) * dx).array();
	Array3d hp = wp - b.p(ip);
	if (!hp.unaryExpr([](double h) { return is_wet(h); }).all())
		gradient.setZero();
	return gradient;
}

void MUSCLObject::full_wet_reconst(index i) {
	assert(full_wet_cell(i));
	auto m = b.mesh();
	auto const& ip = m.triang_points(i);
	auto const& ie = m.triang_edges(i);
	auto const& it = m.triang_triangs(i);
	auto df = gradient(i);
	dw_.col(i) = df.row(0).transpose().array();
	for (int k = 0; k < 3; k++) {
		Edg.prim(ie[k], i, it[k]) = Vol.prim(i) + (df * (m.e(ie[k]) - m.t(i)).matrix()).array();
		double wpk = Vol.w(i) + df.row(0) * (m.p(ip[k]) - m.t(i)).matrix();
		max_wp_[ip[k]] = std::max(max_wp_[ip[k]], wpk);
	}
}