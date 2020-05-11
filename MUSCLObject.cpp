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
} // anonym namespace

namespace utils {
	MUSCLObject::MUSCLObject(VolumeField&& v0)
		: Vol_(std::move(v0))
		, Edg_(Vol_.bathymetry(), 2 * Vol_.bathymetry().mesh().num_edges())
	{}
	bool MUSCLObject::is_dry_cell(index i) const {
		return !is_wet(Vol_.h(i));
	}
	bool MUSCLObject::is_full_wet_cell(index i) const {
		auto const& m = bathymetry().mesh();
		double max_bp = m.triang_points(i).unaryExpr([this](index i) {
			return bathymetry().p(i);
		}).maxCoeff();
		return max_bp >= Vol_.w(i) ? false : m.triang_triangs(i).unaryExpr([this](index i) {
			return i >= 0 ? !is_dry_cell(i) : false;
		}).all();
	}
	bool MUSCLObject::is_part_wet_cell(index i) const {
		return !is_dry_cell(i) && !is_full_wet_cell(i);
	}

	Matrix32d MUSCLObject::gradient(index i) const {
		Matrix32d gradient(Matrix32d::Zero());
		if (!is_full_wet_cell(i))
			return gradient;
		auto const& m = bathymetry().mesh();
		auto const& it = m.triang_triangs(i);
		auto const& ie = m.triang_edges(i);
		Matrix32d points; // points needed to build reconstruction plane 
		Matrix3d  values; // values at these points
		for (int k = 0; k < 3; k++) {
			if (is_full_wet_cell(it[k])) {
				points.row(k) = m.t(it[k]);
				values.col(k) = Vol_.prim(it[k]).matrix();
			}
			else {
				//points.row(k) = m.e(ie[k]);
				//values.col(k) = 0.5 * (Vol_.prim(i) + Edg_.prim(ie[k], it[k], i)).matrix();
				return gradient;
			}
		}
		gradient = values * gradient_coefs(points);
		auto const& ip = m.triang_points(i);
		auto dx = (m.p(ip) - m.t(i).replicate<1, 3>()).matrix();
		Array3d wp = Vol_.w(i) + (gradient.row(0) * dx).array();
		Array3d hp = wp - bathymetry().p(ip);
		if (!hp.unaryExpr([](double h) { return is_wet(h); }).all())
			gradient.setZero();
		return gradient;
	}

	void MUSCLObject::reconstruct_full_wet_cell(index i) {
		auto const& m = bathymetry().mesh();
		auto const& ip = m.triang_points(i);
		auto const& ie = m.triang_edges(i);
		auto const& it = m.triang_triangs(i);
		auto df = gradient(i);
		for (int k = 0; k < 3; k++) {
			Edg_.prim(ie[k], i, it[k]) = Vol_.prim(i) + (df * (m.e(ie[k]) - m.t(i)).matrix()).array();
			double wpk = Vol_.w(i) + df.row(0) * (m.p(ip[k]) - m.t(i)).matrix();
			Max_W_[ip[k]] = std::max(Max_W_[ip[k]], wpk);
		}
	}

	void MUSCLObject::reconstruct_all() {
		Max_W_ = bathymetry().buffer();
		auto const& m = bathymetry().mesh();
		for (size_t i = 0; i < m.num_triangles(); ++i)
			if (is_full_wet_cell(i))
				reconstruct_full_wet_cell(i);
	}
} // namespace utils