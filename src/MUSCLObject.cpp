#include "Includes.h"
#include <CubicPolyMath.h>
#include <MUSCLObject.h>
#include <PointOperations.h>
#include <fstream>
#include <iostream>

MUSCLObject::MUSCLObject(const Domain& b, const VolumeField& v0)
	: m_b(b)
	, m_vol(b, v0)
{}

bool MUSCLObject::IsDryCell(Idx i) const {
	return !IsWet(m_vol.h(i));
}

bool MUSCLObject::IsFullWetCell(Idx i) const {
	const auto& b = GetDomain();
    const auto& m = b.GetTopology();
	if (m.IsTriangleBoundary(i)) {
		return false;
    }
	double max_bp = b.P(m.TriangPoints(i)).row(2).maxCoeff();
	return max_bp < m_vol.w(i);
}

bool MUSCLObject::IsPartWetCell(Idx i) const {
	return !IsDryCell(i) && !IsFullWetCell(i);
}

MUSCLObject::MUSCL MUSCLObject::ReconstructDryCell(Idx i) const {
 	const auto& b = GetDomain();
	Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
    grad.row(0) = b.TriangSlope(i).transpose();
    return MUSCL{b, Array<3>{m_vol.b(i), 0., 0.}, grad, i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructFullWetCell(Idx i) const {
	const auto& b = GetDomain();
	const auto& m = b.GetTopology();
	const auto& ip = m.TriangPoints(i);
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i); 

	Eigen::Matrix3d grad_points; // points needed to build reconstruction plane 
	Eigen::Matrix3d grad_values; // values at these points
	
	for (short int k = 0; k < 3; ++k) {
		if (IsFullWetCell(it[k])) {
	        grad_points.col(k) = b.T(it[k]);
			grad_values.col(k) = m_vol.prim(it[k]).matrix();
		} else if (IsDryCell(it[k])) {
            return MUSCL{b, m_vol.prim(i), Eigen::Matrix32d::Zero(), i};
        } else {
            const auto muscl = ReconstructPartWetCell1(it[k]);
            const auto pt = b.E(ie[k]);
			grad_points.col(k) = pt;
			grad_values.col(k) = 0.5 * (m_vol.prim(i) + muscl.AtPoint(pt)).matrix();
		}
	}


	Eigen::Matrix32d df = Eigen::Matrix32d::Zero();
    df.row(0) = Gradient(grad_points);

	const auto& dx = (b.P(ip).matrix() * (Eigen::Matrix3d::Identity() - Eigen::Matrix3d::Constant(1. / 3.))).topRows(2);
	Array<3> wp = (df.row(0) * dx).array() + m_vol.w(i);
	Array<3> hp = wp - b.P(ip);

	if (!hp.unaryExpr([](double h) { return IsWet(h); }).all()) {
		df.setZero();
    }

	Array<3> TVD = Array<3>::Constant(1.);

	for (int k = 0; k < 3; k++) {
		Array<3> vtmin = m_vol.prim(i).min(m_vol.prim(it[k]));
		Array<3> vtmax = m_vol.prim(i).max(m_vol.prim(it[k]));	
		Array<3> vek   = m_vol.prim(i) + (df * (b.E(ie[k]) - b.T(i)).matrix().topRows(2)).array();
		TVD = ((vtmin <= vek) && (vek <= vtmax)).select(TVD, Array<3>::Zero());
	}
	
    return MUSCL{b, m_vol.prim(i), TVD.matrix().asDiagonal() * df, i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructPartWetCell1(Idx i) const {
    const auto& b = GetDomain();
	const auto& m = b.GetTopology();
	const auto& ip = m.TriangPoints(i);
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i);

	double b13 = b.P(ip).row(2).maxCoeff();
	double b23 = b.P(ip).row(2).minCoeff();
	double b12 = 3. * m_vol.b(i) - b23 - b13;

	double b_delimiter = b12 + (1./3.) * (b13 - b12) * (b13 - b12) / (b13 - b23);
	double w_rec;
	
    if (m_vol.w(i) >= b13) {
        w_rec = m_vol.w(i);
    } else if (m_vol.w(i) <= b_delimiter) {
		w_rec = b23 + cbrt(3. * m_vol.h(i) * (b13 - b23) * (b12 - b23));
	} else {
		double a = -3.*b13;
		double b =  3.*(b12 * b13 + b13 * b23 - b12 * b23);
		double c = (b13 - b23) * (3. * m_vol.h(i) * (b13 - b12) - b12 * (b12 + b23)) - b23 * b23 * b13;
		w_rec = Bisection(CubicPoly(c, b, a), /* find from */ b12, /* to */ b13);
	}

    return MUSCL{b, Array<3>{w_rec, m_vol.u(i), m_vol.v(i)}, Eigen::Matrix32d::Zero(), i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructPartWetCell2(Idx i) const {
    const auto& b = GetDomain();
  	const auto& m = b.GetTopology();
	
    auto ip = m.TriangPoints(i);
    if (b.P(ip[0])[2] > b.P(ip[1])[2]) std::swap(ip[0], ip[1]);
    if (b.P(ip[1])[2] > b.P(ip[2])[2]) std::swap(ip[1], ip[2]);
    if (b.P(ip[0])[2] > b.P(ip[1])[2]) std::swap(ip[0], ip[1]);

    double b23 = b.P(ip[0])[2];
	double b12 = b.P(ip[1])[2];
    double b13 = b.P(ip[2])[2];

    if ((m_vol.w(i) > b13) || (b13 - b23 < tol)) {
        return ReconstructPartWetCell1(i);
    }

    double w23 = m_max_wp[ip[0]];
    double h23 = w23 - b23;
    double ratio_b = (b12 - b23) / (b13 - b23);
  
    double h_delimiter1 = 1./3. * h23 * ratio_b;
    double h_delimiter2 = 1./3. * h23 * (2. * b13 - b12 - b23) / (b13 - b23);

    double hi = m_vol.h(i);
    double ratio_h = hi / h23;

    Eigen::Matrix3d points;
    points.col(0) = b.P(ip[0]);
    points(0, 2) = w23;

    if (hi <= h_delimiter1) { // 1 point wet, 2 points dry
        
        double k2 = sqrt(3. * ratio_h / ratio_b);
        points.col(1) = k2 * b.P(ip[1]) + (1. - k2) * b.P(ip[0]);
        
        double k3 = sqrt(3. * ratio_h * ratio_b);
        points.col(2) = k3 * b.P(ip[2]) + (1. - k3) * b.P(ip[0]);
    
    } else if (hi >= h_delimiter2) { // 3 points wet
        
        double delta_w = 1.5 * (hi - h_delimiter2);
        
        points.col(1) = b.P(ip[1]);
        points(1, 2) += delta_w;
        points(1, 2) += (1. - ratio_b) * h23;

        points.col(2) = b.P(ip[2]);
        points(2, 2) += delta_w;
    
    } else { // 2 points wet, 1 point dry
        
        double alpha = 3. * ratio_h;
        double beta  = (b13 - b12) / (b13 - b23);
        
        double k1 = 1. - Bisection(CubicPoly{(1. + beta - alpha)/(beta*beta), (alpha - 3.)/beta});
        if (k1 < tol) {
            return ReconstructPartWetCell1(i);
        }
        
        double k3 = 1. - beta * (1. - k1);
        
        double bp1 = k1 * b13 + (1. - k1) * b12;
        double bp3 = k3 * b13 + (1. - k3) * b23;

        points.col(1) = b.P(ip[1]);
        points(1, 2) = b12 + (k1 / k3) * beta * h23;
 
        points.col(2) = k1 * b.P(ip[2]) + (1. - k1) * b.P(ip[1]);
        points(2, 2) = bp1;
    }
 
    Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
    grad.row(0) = Gradient(points).transpose();

    double wt = w23 + grad.row(0) * (b.T(i) - b.P(ip[0])).topRows(2).matrix();
    return MUSCL{b, Array<3>{wt, m_vol.u(i), m_vol.v(i)}, grad, i};
}
