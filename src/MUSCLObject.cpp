#include <PointOperations.h>
#include <MUSCLObject.h>
#include <CubicPolyMath.h>
#include <Solver.h>
#include <iostream>
#include <fstream>

MUSCLObject::MUSCLObject(Bathymetry&& b, const VolumeField& v0)
	: b_(std::move(b))
	, Vol_(b_, v0)
	, Edg_(b_, 2 * b_.mesh().num_edges())
	, S_  (b_, b_.mesh().num_triangles())
{
	std::cout << Vol_.storage_size() << '\n';
}

MUSCLObject::MUSCLObject(Bathymetry&& b, VolumeField&& v0)
	: b_(std::move(b))
	, Vol_(b_, std::move(v0))
	, Edg_(b_, 2 * b_.mesh().num_edges())
	, S_  (b_, 2 * b_.mesh().num_triangles())
{
	std::cout << Vol_.storage_size() << '\n';
}

void MUSCLObject::Dump(const std::string& filename) const {
	std::ofstream fout(filename);
	fout << "x\ty\th\tw\tdwx\tdwy\tu\tv\thu\thv\tb\tis full. wet?\tis part. wet?\tis dry?\n";
	
	for (size_t i = 0; i < mesh().num_triangles(); i++) {
		
		const auto& t = mesh().t(i);
		const auto& ie = mesh().triang_edges(i);
		const auto& it = mesh().triang_triangs(i);

		fout << t[0] << '\t' << t[1] << '\t';
		fout << Vol_.h(i) << '\t' << Vol_.w(i) << '\t';

		fout << S_.u(i) << '\t' << S_.v(i) << '\t';

		fout << Vol_.u(i) << '\t' << Vol_.v(i) << '\t'; 
		fout << Vol_.hu(i) << '\t' << Vol_.hv(i) << '\t';
		fout << b_.at_cell(i) << '\t';

		if (is_full_wet_cell(i))
			fout << "1\t0\t0";
		else if (is_dry_cell(i))
			fout << "0\t0\t1";
		else
			fout << "0\t1\t0";

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

Eigen::Matrix32d MUSCLObject::gradient(idx i) const {

	Eigen::Matrix32d gradient(Eigen::Matrix32d::Zero());

	if (!is_full_wet_cell(i))
		return gradient;

	auto const& m = mesh();
	auto const& it = m.triang_triangs(i);
	auto const& ie = m.triang_edges(i);

	Eigen::Matrix32d points; // points needed to build reconstruction plane 
	Eigen::Matrix3d  values; // values at these points
	
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
	Eigen::Matrix23d dx = (m.p(ip) - m.t(i).replicate<1, 3>()).matrix();
	Array<3> wp = (gradient.row(0) * dx).array() + Vol_.w(i);
	Array<3> hp = wp - bath().at_nodes(ip).transpose();
	if (!hp.unaryExpr([](double h) { return is_wet(h); }).all())
		gradient.setZero();

	return gradient;
}

void MUSCLObject::reconstruct_dry_cell(idx i) {

	auto const& m = mesh();
	auto const& ie = m.triang_edges(i);
	auto const& it = m.triang_triangs(i);

	for (int k = 0; k < 3; k++) {
		Edg_.cons(ie[k], i, it[k]) = Array<3>{Vol_.h(i), 0., 0.};
	}

	S_.vel(i) = Array<2>{0., 0.};
}

void MUSCLObject::reconstruct_full_wet_cell(idx i) {

	auto const& m = mesh();
	auto const& ip = m.triang_points(i);
	auto const& ie = m.triang_edges(i);
	auto const& it = m.triang_triangs(i);

	auto const& Vol = Vol_;

	Eigen::Matrix32d df = gradient(i);
	Array<3> TVD = Array<3>::Constant(1.);

	for (int k = 0; k < 3; k++) {
		Array<3> Vtmin = Vol.prim(i).min(Vol.prim(it[k]));
		Array<3> Vtmax = Vol.prim(i).max(Vol.prim(it[k]));	
		Array<3> Vek = Vol.prim(i) + (df * (m.e(ie[k]) - m.t(i)).matrix()).array();
		TVD = ((Vtmin <= Vek) && (Vek <= Vtmax)).select(TVD, Array<3>::Zero());
	}
	
	Eigen::Matrix32d df_lim = TVD.matrix().asDiagonal() * df;

	for (int k = 0; k < 3; k++) {
		Array<3> Vek = Vol.prim(i) + (df_lim * (m.e(ie[k]) - m.t(i)).matrix()).array();
		Edg_.prim(ie[k], i, it[k]) = Vek;
		double wpk = Vol.w(i) + df_lim.row(0) * (m.p(ip[k]) - m.t(i)).matrix();
		Max_W_[ip[k]] = std::max(Max_W_[ip[k]], wpk);
	}

	S_.vel(i) = df_lim.row(0).transpose().array();
}

void MUSCLObject::reconstruct_part_wet_cell(idx i) {

	auto const& m = mesh();
	auto const& b = bath();
	auto const& ip = m.triang_points(i);
	auto const& ie = m.triang_edges(i);
	auto const& it = m.triang_triangs(i);

	double B13 = b.at_nodes(ip).maxCoeff();
	double B23 = b.at_nodes(ip).minCoeff();
	double B12 = 3. * b.at_cell(i) - B23 - B13;

	double B_delimiter = B12 + (1./3.) * (B13 - B12) * (B13 - B12) / (B13 - B23);
	double w_rec;
	
	if (Vol_.w(i) <= B_delimiter) {
		w_rec = B23 + cbrt(3. * Vol_.h(i) * (B13 - B23) * (B12 - B23));
	} else {
		double coef_2nd_degree = -3.*B13;
		double coef_1st_degree = 3.*(B12 * B13 + B13 * B23 - B12 * B23);
		double coef_free = (B13 - B23) * (3. * Vol_.h(i) * (B13 - B12) - B12 * (B12 + B23))
										 - B23 * B23 * B13;
		w_rec = bisection(
			cubic::poly(coef_2nd_degree, coef_1st_degree, coef_free),
			/* find from */ B12, /* to */ B13);
	}

	for (int k = 0; k < 3; k++) {
		Max_W_[ip[k]] = std::max(Max_W_[ip[k]], w_rec);
		double shore = b.at_edge(ie[k]);
		if (w_rec > shore) {
			Edg_.prim(ie[k], i, it[k]) = Array<3>{w_rec, Vol_.u(i), Vol_.v(i)};
		} else {
			Edg_.prim(ie[k], i, it[k]) = Array<3>{shore, 0, 0};
		}
	}

	S_.vel(i) = Array<2>{0., 0.};
}


void MUSCLObject::reconstruct_all() {
		
	Max_W_ = bath().buff();
	auto const& m = mesh();
	for (size_t i = 0; i < m.num_triangles(); ++i) {
		if (is_full_wet_cell(i))
			reconstruct_full_wet_cell(i);
		else if (is_dry_cell(i))
			reconstruct_dry_cell(i);
		else
			reconstruct_part_wet_cell(i);
	}
}
