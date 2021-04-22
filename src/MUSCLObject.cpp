#include <CubicPolyMath.h>
#include <MUSCLObject.h>
#include <PointOperations.h>
#include <fstream>
#include <iostream>

MUSCLObject::MUSCLObject(Bathymetry&& b, const VolumeField& v0)
	: m_b(std::move(b))
	, m_vol(m_b, v0)
	, m_edg(m_b, 2 * m_b.Mesh().NumEdges())
	//, m_s  (m_b,     m_b.Mesh().NumTriangles())
  , m_src(m_b, 2 * m_b.Mesh().NumEdges())
{}

MUSCLObject::MUSCLObject(Bathymetry&& b, VolumeField&& v0)
	: m_b(std::move(b))
	, m_vol(m_b, std::move(v0))
	, m_edg(m_b, 2 * m_b.Mesh().NumEdges())
	//, m_s  (m_b, 2 * m_b.Mesh().NumTriangles())
  , m_src(m_b, 2 * m_b.Mesh().NumEdges())
{}

void MUSCLObject::DumpFields(const std::string& filename) const {

	std::ofstream fout(filename);
	fout << "x\ty\th\tw\tdwx\tdwy\tu\tv\thu\thv\tb\tis full. wet?\tis part. wet?\tis dry?\n";
	
	for (size_t i = 0; i < Mesh().NumTriangles(); i++) {
		
		const auto& t = Mesh().T(i);
		const auto& ie = Mesh().TriangEdges(i);
		const auto& it = Mesh().TriangTriangs(i);

		fout << t[0] << '\t' << t[1] << '\t';

		fout << m_vol.h(i) << '\t' << m_vol.w(i) << '\t';

		//fout << m_s.u(i) << '\t' << m_s.v(i) << '\t';

		fout << m_vol.u(i)  << '\t' << m_vol.v(i) << '\t'; 
		fout << m_vol.hu(i) << '\t' << m_vol.hv(i) << '\t';
		fout << m_vol.b(i)  << '\t';

		if (IsFullWetCell(i))
			fout << "1\t0\t0";
		else if (IsDryCell(i))
			fout << "0\t0\t1";
		else
			fout << "0\t1\t0";

		fout << '\n';
	}
	fout.close();
}

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

Eigen::Matrix32d MUSCLObject::Gradient(Idx i) const {

  Eigen::Matrix32d gradient(Eigen::Matrix32d::Zero());

	if (!IsFullWetCell(i))
		return gradient;

	const auto& m = Mesh();
	const auto& it = m.TriangTriangs(i);
	const auto& ie = m.TriangEdges(i);

	Eigen::Matrix32d points; // points needed to build reconstruction plane 
	Eigen::Matrix3d  values; // values at these points
	
	for (int k = 0; k < 3; k++) {
		if (IsFullWetCell(it[k])) {
			points.row(k) = m.T(it[k]);
			values.col(k) = m_vol.prim(it[k]).matrix();
		}
		else if (IsDryCell(it[k])) {
      return gradient;
    }
    else {
      
			points.row(k) = m.E(ie[k]);
			values.col(k) = 0.5 * (m_vol.prim(i) + m_edg.prim(ie[k], it[k], i)).matrix();
      
      //return gradient;
		}
	}

	gradient = values * GradientCoefs(points);

	const auto& ip = m.TriangPoints(i);
	Eigen::Matrix23d dx = (m.P(ip) - m.T(i).replicate<1, 3>()).matrix();

	Array<3> wp = (gradient.row(0) * dx).array() + m_vol.w(i);
	Array<3> hp = wp - Bath().AtNodes(ip).transpose();

	if (!hp.unaryExpr([](double h) { return IsWet(h); }).all())
		gradient.setZero();

	return gradient;
}

void MUSCLObject::ReconstructDryCell(Idx i) {

	const auto& m = Mesh();
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i);

	for (int k = 0; k < 3; k++) {
		m_edg.prim(ie[k], i, it[k]) = DryState(m_edg.b(ie[k], i, it[k]));
    m_src.vel(ie[k], i, it[k]) = Bath().Gradient(i);
  }

  //m_s.vel(i) = Array<2>{0.,0.};
}

void MUSCLObject::ReconstructFullWetCell(Idx i) {

	const auto& m = Mesh();
	const auto& ip = m.TriangPoints(i);
	const auto& ie = m.TriangEdges(i);
	const auto& it = m.TriangTriangs(i);

	const auto& vol = m_vol;

	Eigen::Matrix32d df = Gradient(i);
	Array<3> TVD = Array<3>::Constant(1.);

	for (int k = 0; k < 3; k++) {
		Array<3> vtmin = vol.prim(i).min(vol.prim(it[k]));
		Array<3> vtmax = vol.prim(i).max(vol.prim(it[k]));	
		Array<3> vek   = vol.prim(i) + (df * (m.E(ie[k]) - m.T(i)).matrix()).array();
		TVD = ((vtmin <= vek) && (vek <= vtmax)).select(TVD, Array<3>::Zero());
	}
	
	Eigen::Matrix32d df_lim = TVD.matrix().asDiagonal() * df;

	for (int k = 0; k < 3; k++) {
		
    double wpk = vol.w(i) + df_lim.row(0) * (m.P(ip[k]) - m.T(i)).matrix();
    m_max_w[ip[k]] = std::max(m_max_w[ip[k]], wpk);

    Array<3> vek = vol.prim(i) + (df_lim * (m.E(ie[k]) - m.T(i)).matrix()).array();
    m_edg.prim(ie[k], i, it[k]) = vek;
		m_src.vel(ie[k], i, it[k]) = df_lim.row(0).transpose().array();
	}

	//m_s.vel(i) = df_lim.row(0).transpose().array();
}

void MUSCLObject::ReconstructPartWetCell1(Idx i) {

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
	
	if (m_vol.w(i) <= b_delimiter) {

		w_rec = b23 + cbrt(3. * m_vol.h(i) * (b13 - b23) * (b12 - b23));

	} else {

		double a = -3.*b13;
		double b =  3.*(b12 * b13 + b13 * b23 - b12 * b23);
		double c = (b13 - b23) * (3. * m_vol.h(i) * (b13 - b12) - b12 * (b12 + b23)) - b23 * b23 * b13;

		w_rec = Bisection(CubicPoly(c, b, a), /* find from */ b12, /* to */ b13);
	}

	for (int k = 0; k < 3; k++) {
		
    m_max_w[ip[k]] = std::max(m_max_w[ip[k]], w_rec);
		
    double shore = m_edg.b(ie[k], i, it[k]);
		if (w_rec > shore) {
			m_edg.prim(ie[k], i, it[k]) = Array<3>{w_rec, m_vol.u(i), m_vol.v(i)};
		} else {
			m_edg.prim(ie[k], i, it[k]) = Array<3>{shore, 0, 0};
		}
	}
}

void MUSCLObject::ReconstructPartWetCell2(Idx i) {
  
  const auto& m = Mesh();
	const auto& b = Bath();

  // We sort points:
  // 1) "lowest" point ip[0] 
  // 2) "middle" point ie[1]
  // 3) "highest" point ip[2]

  auto ip = m.TriangPoints(i);
  
  if (b.AtNode(ip[0]) > b.AtNode(ip[1])) {
    std::swap(ip[0], ip[1]);
  }
  if (b.AtNode(ip[1]) > b.AtNode(ip[2])) {
    std::swap(ip[1], ip[2]);
  }
  if (b.AtNode(ip[0]) > b.AtNode(ip[1])) {
    std::swap(ip[0], ip[1]);
  }

  double b23 = b.AtNode(ip[0]);
	double b12 = b.AtNode(ip[1]);
  double b13 = b.AtNode(ip[2]);

  double w23 = m_max_w[ip[0]];
  double h23 = w23 - b23;
  double ratio_b = (b12 - b23) / (b13 - b23);
  
  double h_delimiter1 = 1./3. * h23 * ratio_b;
  double h_delimiter2 = 1./3. * h23 * (2. - ratio_b);

  double hi = m_vol.h(i);
  double ratio_h = hi / h23;

  //Further, we will distinguish 3 kinds of edges:
  // 1) "lowest" edge ie[0] 
  // 2) "middle" edge ie[1]
  // 3) "highest" one ie[2]

  auto ie = m.TriangEdges(i);
  auto it = m.TriangTriangs(i);

  if (m_edg.b(ie[0], i, it[0]) > m_edg.b(ie[1], i, it[1])) {
    std::swap(ie[0], ie[1]);
    std::swap(it[0], it[1]);
  }
  if (m_edg.b(ie[1], i, it[1]) > m_edg.b(ie[2], i, it[2])) {
    std::swap(ie[1], ie[2]);
    std::swap(it[1], it[2]);
  }
  if (m_edg.b(ie[0], i, it[0]) > m_edg.b(ie[1], i, it[1])) {
    std::swap(ie[0], ie[1]);
    std::swap(it[0], it[1]);
  }

  Eigen::Matrix13d grad_values;
  Eigen::Matrix32d grad_points;
  
  grad_values[0] = w23;
  grad_points.row(0) = m.P(ip[0]);

  auto semiwet_edge = [this, &i, &ie, &it]
    (int k, double val1, double val2, double portion)
  {
    //if (portion < 0.5) {
    //  return DryState(m_edg.b(ie[k], i, it[k]));
    //}
    //else {
      double c = 0.5 / portion;
      double w_rec = c * val1 + (1. - c) * val2;
      return Array<3>{ w_rec, m_vol.u(i), m_vol.v(i) };
    //}
  };

  //Case 1
  if (hi <= h_delimiter1) {

    m_edg.prim(ie[2], i, it[2]) = DryState(m_edg.b(ie[2], i, it[2]));
    m_src.vel(ie[2], i, it[2])  = Bath().Gradient(i);

    double k2 = sqrt(3. * ratio_h / ratio_b);
    grad_values    [1] = k2 * b.AtNode(ip[1]) + (1. - k2) * b.AtNode(ip[0]);
    grad_points.row(1) = k2 * m.P(ip[1]) + (1. - k2) * m.P(ip[0]);
    
    double k3 = sqrt(3. * ratio_h * ratio_b);
    grad_values    [2] = k3 * b.AtNode(ip[2]) + (1. - k3) * b.AtNode(ip[0]);
    grad_points.row(2) = k3 * m.P(ip[2]) + (1. - k3) * m.P(ip[0]);
      
    Array<2> gradient = (grad_values * GradientCoefs(grad_points)).transpose().array();

    if (k2 < 0.5) {
      m_edg.prim(ie[0], i, it[0]) = DryState(m_edg.b(ie[0], i, it[0]));
      m_src.vel(ie[0], i, it[0])  = Bath().Gradient(i);
    } else {
      m_edg.prim(ie[0], i, it[0]) = semiwet_edge(0, grad_values[0], grad_values[1], k2);
      m_src.vel(ie[0], i, it[0]) = gradient;
    }

    if (k3 < 0.5) {
      m_edg.prim(ie[1], i, it[1]) = DryState(m_edg.b(ie[1], i, it[1]));
      m_src.vel(ie[1], i, it[1])  = Bath().Gradient(i);
    } else {
      m_edg.prim(ie[1], i, it[1]) = semiwet_edge(1, grad_values[0], grad_values[2], k3);
      m_src.vel(ie[1], i, it[1]) = gradient;
    }
 }
  //Case 3
  else if (hi >= h_delimiter2) {

    double delta_w = 1.5 * (hi - h_delimiter2);

    grad_values[1] = delta_w + b.AtNode(ip[1]) + (1. - ratio_b) * h23;
    grad_points.row(1) = m.P(ip[1]);

    grad_values[2] = delta_w + b.AtNode(ip[2]);
    grad_points.row(2) = m.P(ip[2]);

    double w_rec;

    w_rec = 0.5 * (grad_values[0] + grad_values[1]);
    m_edg.prim(ie[0], i, it[0]) = Array<3>{w_rec, m_vol.u(i), m_vol.v(i)};
    w_rec = 0.5 * (grad_values[0] + grad_values[2]);
    m_edg.prim(ie[1], i, it[1]) = Array<3>{w_rec, m_vol.u(i), m_vol.v(i)};
    w_rec = 0.5 * (grad_values[1] + grad_values[2]);
    m_edg.prim(ie[2], i, it[2]) = Array<3>{w_rec, m_vol.u(i), m_vol.v(i)};

    Array<2> gradient = (grad_values * GradientCoefs(grad_points)).transpose().array();
    for (short int k = 0; k < 3; k++)
      m_src.vel(ie[k], i, it[k]) = gradient;
  }
  //Case 2
  else {

    double alpha = 3. * ratio_h;
    double beta  = 1. - ratio_b;

    double k1 = 1. - Bisection(CubicPoly{(1. + beta - alpha)/(beta*beta), (alpha - 3.)/beta});
    double k3 = 1. - beta * (1. - k1);

    double bp1 = k1 * b.AtNode(ip[2]) + (1. - k1) * b.AtNode(ip[1]);
    double bp3 = k3 * b.AtNode(ip[2]) + (1. - k3) * b.AtNode(ip[0]);

    grad_values    [1] = b.AtNode(ip[1]) + (k1 / k3) * beta * h23;
    grad_points.row(1) = m.P(ip[1]);

    grad_values    [2] = bp1;
    grad_points.row(2) = k1 * m.P(ip[2]) + (1. - k1) * m.P(ip[1]);
    
    Array<2> gradient = (grad_values * GradientCoefs(grad_points)).transpose().array();
    double w_rec = 0.5 * (grad_values[0] + grad_values[1]);

    m_edg.prim(ie[0], i, it[0]) = Array<3>{w_rec, m_vol.u(i), m_vol.v(i)};
    m_src.vel(ie[0], i, it[0]) = gradient;
  
    if (k3 < 0.5) {
      m_edg.prim(ie[1], i, it[1]) = DryState(m_edg.b(ie[1], i, it[1]));
      m_src.vel(ie[1], i, it[1]) = Bath().Gradient(i);
    } else {
      m_edg.prim(ie[1], i, it[1]) = semiwet_edge(1, grad_values[0], bp3, k3);
      m_src.vel(ie[1], i, it[1]) = gradient;
    }

    if (k1 < 0.5) {
      m_edg.prim(ie[2], i, it[2]) = DryState(m_edg.b(ie[2], i, it[2]));
      m_src.vel(ie[2], i, it[2]) = Bath().Gradient(i);
    } else {
      m_edg.prim(ie[2], i, it[2]) = semiwet_edge(2, grad_values[1], bp1, k1);
      m_src.vel(ie[2], i, it[2]) = gradient;
    }
  }
 
  //m_s.vel(i) = (grad_values * GradientCoefs(grad_points)).transpose().array();
}
	
void MUSCLObject::ReconstructAll() {
		
	m_max_w = Bath().Buff();
	
  const auto& m = Mesh();
	
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		if (IsDryCell(i))
			ReconstructDryCell(i);
		else if (!IsFullWetCell(i))
			ReconstructPartWetCell1(i);
	}

  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		if (IsFullWetCell(i))
			ReconstructFullWetCell(i);
	}

  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		if (IsPartWetCell(i))
			ReconstructPartWetCell2(i);
	}

}
