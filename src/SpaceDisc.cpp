#include <SpaceDisc.h>

SpaceDisc::SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, const VolumeField& v0, double cor, double tau)
  : MUSCLObject(std::move(b), v0)
  , m_edg(m_b, 2 * m_b.Mesh().NumEdges())
  , m_src(m_b, 2 * m_b.Mesh().NumEdges())
  , m_fluxer(fluxer)
  , m_cor(cor)
  , m_tau(tau)
{
	m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
}

SpaceDisc::SpaceDisc(const Fluxer& fluxer, Bathymetry&& b, VolumeField&& v0, double cor, double tau)
  : MUSCLObject(std::move(b), std::move(v0))
  , m_edg(m_b, 2 * m_b.Mesh().NumEdges())
  , m_src(m_b, 2 * m_b.Mesh().NumEdges())
  , m_fluxer(fluxer)
  , m_cor(cor)
  , m_tau(tau)
{
	m_f.resize(Eigen::NoChange, this->Mesh().NumEdges());
}

void SpaceDisc::UpdateInterfaceValues(const MUSCL& muscl) {

  const Idx i = muscl.OriginId();
  const auto& m = Mesh();
  const auto& ip = m.TriangPoints(i);
  const auto& ie = m.TriangEdges(i);
  const auto& it = m.TriangTriangs(i);

  for (short int k = 0; k < 3; ++k) {
    m_max_wp[ip[k]] = std::max(m_max_wp[ip[k]], muscl.AtPoint(m.P(ip[k]))[0]);
    m_edg.prim(ie[k], i, it[k]) = muscl.AtPoint(m.E(ie[k]));
    m_src.vel(ie[k], i, it[k]) = muscl.Gradient(m.E(ie[k])).row(0).transpose().array();
    double u = m_edg.u(ie[k], i, it[k]);
    double v = m_edg.v(ie[k], i, it[k]);
    m_src.vel(ie[k], i, it[k]) += m_cor * Array<2>{ -v, u };
  }
}

void SpaceDisc::ComputeInterfaceValues() {
		
	m_max_wp = Bath().Buff();
  const auto& m = Mesh();
	
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		if (IsDryCell(i))
			UpdateInterfaceValues(ReconstructDryCell(i));
		else if (!IsFullWetCell(i))
			UpdateInterfaceValues(ReconstructPartWetCell1(i));
    else
      UpdateInterfaceValues(ReconstructFullWetCell(i));
	}
 
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		if (IsPartWetCell(i))
			UpdateInterfaceValues(ReconstructPartWetCell2(i));
	}
}

void SpaceDisc::ComputeFluxes() {

  m_min_length_to_wavespeed = 1.;
  const auto& m = Mesh();

  for (size_t i = 0; i < m.NumEdges(); ++i) {

    const auto& it = m.EdgeTriangs(i);
    Idx lf = it[0];
    Idx lt = it[1];

    auto n = m.Norm(i, lf);
    
    switch (lt) {
    case static_cast<Idx>(Boundaries::SOLID_WALL):
      m_f.col(i) = ElemFlux(n, Array<3>{ m_vol.h(lf), 0., 0. });
      break;
    default:
      m_f.col(i) = m_fluxer(this, i, lf, lt, &m_min_length_to_wavespeed);
    }
  }
}

Array<3> SpaceDisc::ComputeIntegrals() {

  Array<3> res = Array<3>::Zero();
  const auto& m = Mesh();
  const size_t n = m.NumTriangles();

  for (size_t i = 0; i < n; ++i) {
    if (IsPartWetCell(i))
      UpdateInterfaceValues(ReconstructPartWetCell1(i));
		if (IsFullWetCell(i)) {
      double h = m_vol.h(i);
      double b = m_vol.b(i);
      double u = m_vol.u(i);
      double v = m_vol.v(i);
      res += m.Area(i) * Array<3>{h, 0.5 * h * (u * u + v * v), 0.5 * h * h + h * b};
			UpdateInterfaceValues(ReconstructFullWetCell(i));
    }
	}

  for (size_t i = 0; i < n; ++i) {
		if (IsPartWetCell(i)) {
      MUSCL muscl = ReconstructPartWetCell2(i);
      const auto& ip = m.TriangPoints(i);
      res[0] += m.Area(i) * m_vol.h(i);
      res.tail<2>() += m.Area(i) * TriangAverage<2,10>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), [&](const Point& p) {
        Array<3> z = muscl.AtPoint(p);
        double b = Bath().AtPoint(muscl.OriginId(), p);
        double h = z[0] - b;
        Array<2> res = Array<2>{ 0.5 * h * (z[1] * z[1] + z[2] * z[2]), 0.5 * h * h + h * b};
        return res;
      });
    }
	}

  return res;
}

Storage<3> SpaceDisc::CompareWith(const Test& test, double t) {

  m_max_wp = Bath().Buff();
  const auto& m = Mesh();
  const size_t n = m.NumTriangles();

  Storage<3> errL2;
  errL2.resize(Eigen::NoChange, n);
  	
  for (size_t i = 0; i < n; ++i) {
		if (IsDryCell(i)) {
      MUSCL muscl = ReconstructDryCell(i);
      const auto& ip = m.TriangPoints(i);
      errL2.col(i) = TriangAverage<3,1>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), [&](const Point& p) {
        double b = Bath().AtPoint(muscl.OriginId(), p);
        Array<3> U = Array<3>{ test.w(p[0], p[1], t), test.hu(p[0], p[1], t), test.hv(p[0], p[1], t) };
        Array<3> V = muscl.AtPoint(p);
        V.tail<2>() *= (V[0] -  b);
        Array<3> res = (U - V).abs();
        return res;
      });
			UpdateInterfaceValues(muscl);
    }
		else if (!IsFullWetCell(i))
			UpdateInterfaceValues(ReconstructPartWetCell1(i));
	}

  for (size_t i = 0; i < n; ++i) {
		if (IsFullWetCell(i)) {
      MUSCL muscl = ReconstructFullWetCell(i);
      const auto& ip = m.TriangPoints(i);
      errL2.col(i) = TriangAverage<3,1>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), [&](const Point& p) {
        double b = Bath().AtPoint(muscl.OriginId(), p);
        Array<3> U = Array<3>{ test.w(p[0], p[1], t), test.hu(p[0], p[1], t), test.hv(p[0], p[1], t) };
        Array<3> V = muscl.AtPoint(p);
        V.tail<2>() *= (V[0] - b);
        Array<3> res = (U - V).abs();
        return res;
      });
			UpdateInterfaceValues(muscl);
    }
	}

  for (size_t i = 0; i < n; ++i) {
		if (IsPartWetCell(i)) {
      MUSCL muscl = ReconstructPartWetCell2(i);
      const auto& ip = m.TriangPoints(i);
      errL2.col(i) = TriangAverage<3,100>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), [&](const Point& p) {
        double b = Bath().AtPoint(muscl.OriginId(), p);
        Array<3> U = Array<3>{ test.w(p[0], p[1], t), test.hu(p[0], p[1], t), test.hv(p[0], p[1], t) };
        Array<3> V = muscl.AtPoint(p);
        V.tail<2>() *= (V[0] - b);
        Array<3> res = (U - V).abs();
        return res;
      });
    }
	}

  return errL2;
}

