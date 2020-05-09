#include "ValueField.h"
#include "PointOperations.h"
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
		if (abs(d) < Solver::tol)
			throw SolverError("Division by zero in gradient coefficient computation");
		res /= d;
		return res;
	}
}

bool ValueField::dry_cell(index i) const {
	return !is_wet(h(i));
}
bool ValueField::full_wet_cell(index i) const {
	auto const& mesh = bathymetry().mesh();
	double max_bp = mesh.triang_points(i).unaryExpr([this](index i) {
		return bathymetry().p(i);
	}).maxCoeff();
	return max_bp < w(i) ? false : mesh.triang_triangs(i).unaryExpr([this](index i) {
		return i >= 0 ? dry_cell(i) : false;
	}).any();
}
bool ValueField::part_wet_cell(index i) const {
	return !dry_cell(i) && !full_wet_cell(i);
}
Matrix32d ValueField::gradient(index i) const {
	assert(full_wet_cell(i));
	auto mesh = bathymetry().mesh();
	triang_tag it = mesh.triang_triangs(i);
	if (!it.unaryExpr([](index i) { return i < 0; }).any())
		return Matrix32d::Zero();
	triang_tag const& ie = mesh.triang_edges(i);
	//first step of reconstruction: 3-stencil plane's gradient 
	Matrix32d points; // points needed to build reconstruction plane 
	Matrix3d  values; // values at these points
	for (int k = 0; k < 3; k++) {
		if (full_wet_cell(it[k])) {
			points.row(k) = mesh.t(it[k]);
			values.col(k) = prim(it[k]).matrix();
		}
		else {
			points.row(k) = mesh.e(ie[k]);
			values.col(k) = (i == mesh.edge_triangs(ie[k])[0]) ?
				0.5 * (prim(i) + Ye_[ie[k]][1]).matrix() :
				0.5 * (prim(i) + Ye_[ie[k]][0]).matrix();
		}
	}
	auto gradient = values * gradient_coefs(points);
	triang_tag const& ip = mesh.triang_points(i);
	auto dx = (mesh.p(ip) - mesh.t(i).replicate<1, 3>()).matrix();
	Array3d wp = w(i) + (gradient.row(0) * dx).array();
	Array3d hp = wp - bathymetry().p(ip);
	if (!hp.unaryExpr([](double h) { return is_wet(h); }).all())
		return Matrix2d::Zero();
	return gradient;
}
