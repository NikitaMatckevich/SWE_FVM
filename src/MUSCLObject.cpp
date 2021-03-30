#include <PointOperations.h>
#include <MUSCLObject.h>
#include <Solver.h>
#include <iostream>
#include <fstream>
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
} // anonymous namespace

MUSCLObject::MUSCLObject(Bathymetry&& b, const VolumeField& v0)
	: b_(std::move(b))
	, Vol_(b_, v0)
	, Edg_(b_)
{
	std::cout << Vol_.storage_size() << '\n';
}

MUSCLObject::MUSCLObject(Bathymetry&& b, VolumeField&& v0)
	: b_(std::move(b))
	, Vol_(b_, std::move(v0))
	, Edg_(b_)
{
	std::cout << Vol_.storage_size() << '\n';
}

void MUSCLObject::Dump(const std::string& filename) const {
	std::ofstream fout(filename);
	fout << "x\ty\tw\tu\tv\thu\thv\tb\n";
	for (size_t i = 0; i < mesh().num_triangles(); i++) {
		const auto& t = mesh().t(i);
		fout << t[0] << '\t' << t[1] << '\t';
		fout << Vol_.w(i) << '\t';
		fout << Vol_.u(i) << '\t' << Vol_.v(i) << '\t'; 
		fout << Vol_.hu(i) << '\t' << Vol_.hv(i) << '\t'; 
		fout << b_.at_cell(i); 
		fout << '\n';
	}
	fout.close();
}

bool MUSCLObject::is_dry_cell(idx i) const {
	return !is_wet(Vol_.h(i));
}

bool MUSCLObject::is_full_wet_cell(idx i) const {
	auto const& m = mesh();
	
	if (m.is_triangle_boundary(i))
		return false;

	double max_bp = bath().at_nodes(m.triang_points(i)).maxCoeff();
	return max_bp < Vol_.w(i);
}

bool MUSCLObject::is_part_wet_cell(idx i) const {
	return !is_dry_cell(i) && !is_full_wet_cell(i);
}

Matrix32d MUSCLObject::gradient(idx i) const {
	Matrix32d gradient(Matrix32d::Zero());
	if (!is_full_wet_cell(i))
		return gradient;
	auto const& m = mesh();
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
	Matrix23d dx = (m.p(ip) - m.t(i).replicate<1, 3>()).matrix();
	Array3d wp = (gradient.row(0) * dx).array() + Vol_.w(i);
	Array3d hp = wp - bath().at_nodes(ip);
	if (!hp.unaryExpr([](double h) { return is_wet(h); }).all())
		gradient.setZero();
	return gradient;
}

void MUSCLObject::reconstruct_full_wet_cell(idx i) {
	auto const& m = mesh();
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
	if (Edg_.storage_size() == 0) {
		Edg_.resize_storage(2 * mesh().num_edges());
	}
	
	Max_W_ = bath().buff();
	auto const& m = mesh();
	for (size_t i = 0; i < m.num_triangles(); ++i) {
		if (is_full_wet_cell(i))
			reconstruct_full_wet_cell(i);
	}
}
