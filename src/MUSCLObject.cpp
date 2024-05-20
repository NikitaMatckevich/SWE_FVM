#include <CubicPolyMath.h>
#include <MUSCLObject.h>
#include <PointOperations.h>
#include <fstream>
#include <iostream>

MUSCLObject::MUSCLObject(Bathymetry&& b, const VolumeField& v0)
	: m_b(std::move(b))
	, m_vol(m_b, v0)
{}

MUSCLObject::MUSCLObject(Bathymetry&& b, VolumeField&& v0)
	: m_b(std::move(b))
	, m_vol(m_b, std::move(v0))
{}

bool MUSCLObject::IsDryCell(Idx i) const {
	return !IsWet(m_vol.h(i));
}

bool MUSCLObject::IsFullWetCell(Idx i) const {
	const auto& m = Mesh();
	if (m.IsTriangleBoundary(i))
		return false;
	double max_bp = Bath().AtNodes(m.TriangPoints(i)).maxCoeff();
	return max_bp < m_vol.w(i);
}

bool MUSCLObject::IsPartWetCell(Idx i) const {
	return !IsDryCell(i) && !IsFullWetCell(i);
}

MUSCLObject::MUSCL MUSCLObject::ReconstructDryCell(Idx i) const {
  Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
  grad.row(0) = Bath().Gradient(i).matrix().transpose();
  return MUSCL{&Bath(), Array<3>{m_vol.b(i), 0., 0.}, grad, i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructFullWetCell(Idx i) const {

	const auto& m = Mesh();
	const auto& ip = m.TriangPoints(i);
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i); 

	Eigen::Matrix32d grad_points; // points needed to build reconstruction plane 
	Eigen::Matrix3d  grad_values; // values at these points
	
	for (short int k = 0; k < 3; ++k) {
		if (IsFullWetCell(it[k])) {
			grad_points.row(k) = m.T(it[k]);
			grad_values.col(k) = m_vol.prim(it[k]).matrix();
		}
		else if (IsDryCell(it[k])) {
      return MUSCL{&Bath(), m_vol.prim(i), Eigen::Matrix32d::Zero(), i};
    }
    else {
      const auto reco = ReconstructPartWetCell1(it[k]);
      const auto pt = m.E(ie[k]);
			grad_points.row(k) = pt;
			grad_values.col(k) = 0.5 * (m_vol.prim(i) + reco.AtPoint(pt)).matrix();
		}
	}

	Eigen::Matrix32d df = grad_values * GradientCoefs(grad_points);
	Eigen::Matrix23d dx = (m.P(ip) - m.T(i).replicate<1, 3>()).matrix();
	Array<3> wp = (df.row(0) * dx).array() + m_vol.w(i);
	Array<3> hp = wp - Bath().AtNodes(ip).transpose();

	if (!hp.unaryExpr([](double h) { return IsWet(h); }).all())
		df.setZero();

	Array<3> TVD = Array<3>::Constant(1.);

	for (int k = 0; k < 3; k++) {
		Array<3> vtmin = m_vol.prim(i).min(m_vol.prim(it[k]));
		Array<3> vtmax = m_vol.prim(i).max(m_vol.prim(it[k]));	
		Array<3> vek   = m_vol.prim(i) + (df * (m.E(ie[k]) - m.T(i)).matrix()).array();
		TVD = ((vtmin <= vek) && (vek <= vtmax)).select(TVD, Array<3>::Zero());
	}
	
  return MUSCL{&Bath(), m_vol.prim(i), TVD.matrix().asDiagonal() * df, i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructPartWetCell1(Idx i) const {

	const auto& m = Mesh();
	const auto& b = Bath();
	const auto& ip = m.TriangPoints(i);
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i);

	double b13 = b.AtNodes(ip).maxCoeff();
	double b23 = b.AtNodes(ip).minCoeff();
	double b12 = 3. * m_vol.b(i) - b23 - b13;

	double b_delimiter = b12 + (1./3.) * (b13 - b12) * (b13 - b12) / (b13 - b23);
	double w_rec;
	
  if (m_vol.w(i) >= b13) {
    w_rec = m_vol.w(i);
  }
	else if (m_vol.w(i) <= b_delimiter) {
		w_rec = b23 + cbrt(3. * m_vol.h(i) * (b13 - b23) * (b12 - b23));
	}
  else {
		double a = -3.*b13;
		double b =  3.*(b12 * b13 + b13 * b23 - b12 * b23);
		double c = (b13 - b23) * (3. * m_vol.h(i) * (b13 - b12) - b12 * (b12 + b23)) - b23 * b23 * b13;
		w_rec = Bisection(CubicPoly(c, b, a), /* find from */ b12, /* to */ b13);
	}

  return MUSCL{&Bath(), Array<3>{w_rec, m_vol.u(i), m_vol.v(i)}, Eigen::Matrix32d::Zero(), i};
}

MUSCLObject::MUSCL MUSCLObject::ReconstructPartWetCell2(Idx i) const {
  
  const auto& m = Mesh();
	const auto& b = Bath();
  
  auto ip = m.TriangPoints(i);
  if (b.AtNode(ip[0]) > b.AtNode(ip[1])) std::swap(ip[0], ip[1]);
  if (b.AtNode(ip[1]) > b.AtNode(ip[2])) std::swap(ip[1], ip[2]);
  if (b.AtNode(ip[0]) > b.AtNode(ip[1])) std::swap(ip[0], ip[1]);

  double b23 = b.AtNode(ip[0]);
	double b12 = b.AtNode(ip[1]);
  double b13 = b.AtNode(ip[2]);

  if ((m_vol.w(i) > b13) || (b13 - b23 < tol))
    return ReconstructPartWetCell1(i);

  double w23 = m_max_wp[ip[0]];
  double h23 = w23 - b23;
  double ratio_b = (b12 - b23) / (b13 - b23);
  
  double h_delimiter1 = 1./3. * h23 * ratio_b;
  double h_delimiter2 = 1./3. * h23 * (2. * b13 - b12 - b23) / (b13 - b23);

  double hi = m_vol.h(i);
  double ratio_h = hi / h23;

  Eigen::Matrix13d grad_values;
  Eigen::Matrix32d grad_points;
  grad_values[0] = w23;
  grad_points.row(0) = m.P(ip[0]);
   
  if (hi <= h_delimiter1) { // Case 1
    double k2 = sqrt(3. * ratio_h / ratio_b);
    grad_values    [1] = k2 * b.AtNode(ip[1]) + (1. - k2) * b.AtNode(ip[0]);
    grad_points.row(1) = k2 * m.P(ip[1]) + (1. - k2) * m.P(ip[0]);
    double k3 = sqrt(3. * ratio_h * ratio_b);
    grad_values    [2] = k3 * b.AtNode(ip[2]) + (1. - k3) * b.AtNode(ip[0]);
    grad_points.row(2) = k3 * m.P(ip[2]) + (1. - k3) * m.P(ip[0]); 
  }
  else if (hi >= h_delimiter2) { // Case 3
    double delta_w = 1.5 * (hi - h_delimiter2);
    grad_values[1] = delta_w + b.AtNode(ip[1]) + (1. - ratio_b) * h23;
    grad_points.row(1) = m.P(ip[1]);
    grad_values[2] = delta_w + b.AtNode(ip[2]);
    grad_points.row(2) = m.P(ip[2]);
  }
  else { // Case 2
    double alpha = 3. * ratio_h;
    double beta  = (b13 - b12) / (b13 - b23);
    double k1 = 1. - Bisection(CubicPoly{(1. + beta - alpha)/(beta*beta), (alpha - 3.)/beta});
    if (k1 < tol) {
      return ReconstructPartWetCell1(i);
    }
    double k3 = 1. - beta * (1. - k1);
    double bp1 = k1 * b.AtNode(ip[2]) + (1. - k1) * b.AtNode(ip[1]);
    double bp3 = k3 * b.AtNode(ip[2]) + (1. - k3) * b.AtNode(ip[0]);
    grad_values    [1] = b.AtNode(ip[1]) + (k1 / k3) * beta * h23;
    grad_points.row(1) = m.P(ip[1]);
    grad_values    [2] = bp1;
    grad_points.row(2) = k1 * m.P(ip[2]) + (1. - k1) * m.P(ip[1]);
  }
 
  Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
  grad.row(0) = grad_values * GradientCoefs(grad_points);
  double wt = w23 + grad.row(0) * (m.T(i) - m.P(ip[0])).matrix();
  return MUSCL{&Bath(), Array<3>{wt, m_vol.u(i), m_vol.v(i)}, grad, i};
}
