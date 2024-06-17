#include <SpaceDisc.h>
#include <memory>

SpaceDisc::SpaceDisc(const Fluxer& fluxer, const Domain& b, const VolumeField& v0, double cor, double tau)
  : MUSCLObject(b, v0)
  , m_edg(m_b, 2 * m_b.GetTopology().NumEdges())
  , m_src(m_b, 2 * m_b.GetTopology().NumEdges())
  , m_fluxer(fluxer)
  , m_cor(cor)
  , m_tau(tau)
{
	m_f.resize(Eigen::NoChange, m_b.GetTopology().NumEdges());
}

void SpaceDisc::UpdateInterfaceValues(const MUSCL& muscl) {
    const Idx i = muscl.OriginId();
    const auto& m = m_b.GetTopology();
    const auto& ip = m.TriangPoints(i);
    const auto& ie = m.TriangEdges(i);
    const auto& it = m.TriangTriangs(i);

    for (short int k = 0; k < 3; ++k) {
        m_max_wp[ip[k]] = std::max(m_max_wp[ip[k]], muscl.AtPoint(m_b.P(ip[k]))[0]);
        const auto& edgeCenter = m_b.E(ie[k]);
        m_edg.prim(ie[k], i, it[k]) = muscl.AtPoint(edgeCenter);
        m_src.vel(ie[k], i, it[k]) = muscl.Gradient(edgeCenter).row(0).transpose().array();
        double u = m_edg.u(ie[k], i, it[k]);
        double v = m_edg.v(ie[k], i, it[k]);
        m_src.vel(ie[k], i, it[k]) += m_cor * Array<2>{ -v, u };
    }
}

void SpaceDisc::ComputeInterfaceValues() {		
    m_max_wp = m_b.GetGeometry().row(2);
    const auto& m = m_b.GetTopology();
	
    for (size_t i = 0; i < m.NumTriangles(); ++i) {
        if (IsDryCell(i)) {
            UpdateInterfaceValues(ReconstructDryCell(i));
        } else if (!IsFullWetCell(i)) {
	        UpdateInterfaceValues(ReconstructPartWetCell1(i));
        } else {
            UpdateInterfaceValues(ReconstructFullWetCell(i));
	    }
    }
 
    for (size_t i = 0; i < m.NumTriangles(); ++i) {
        if (IsPartWetCell(i)) {
	        UpdateInterfaceValues(ReconstructPartWetCell2(i));
        }
    }
}

void SpaceDisc::ComputeFluxes() {

  m_min_length_to_wavespeed = 1.;
  const auto& m = m_b.GetTopology();

  for (size_t i = 0; i < m.NumEdges(); ++i) {
    const auto& it = m.EdgeTriangs(i);
    Idx lf = it[0];
    Idx lt = it[1];

    auto n = m_b.Norm(i, lf);
    
    switch (lt) {
    case static_cast<Idx>(Boundaries::SOLID_WALL):
      m_f.col(i) = ElemFlux(n, Array<3>{ m_vol.h(lf), 0., 0. });
      break;
    default:
      m_f.col(i) = m_fluxer(this, i, lf, lt, &m_min_length_to_wavespeed);
    }
  }
}

/*
Array<3> ComputeIntegrals(const SpaceDisc& sd) {
  Array<3> res = Array<3>::Zero();
  const auto& m = sd.GetTopology();
  const size_t n = m.NumTriangles();

  for (size_t i = 0; i < n; ++i) {
	if (sd.IsFullWetCell(i)) {
      double h = sd.GetVolField().h(i);
      double b = sd.GetVolField().b(i);
      double u = sd.GetVolField().u(i);
      double v = sd.GetVolField().v(i);
      res += m.Area(i) * Array<3>{h, 0.5 * h * (u * u + v * v), 0.5 * h * h + h * b};
    } else if (sd.IsPartWetCell(i)) {
      const MUSCLObject::MUSCL muscl = sd.ReconstructPartWetCell2(i);
      const auto& ip = m.TriangPoints(i);
      res[0] += m.Area(i) * sd.GetVolField().h(i);
      res.tail<2>() += m.Area(i) * TriangAverage<2,10>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), [&](const Point& p) {
        Array<3> z = muscl.AtPoint(p);
        double b = sd.Bath().AtPoint(muscl.OriginId(), p);
        double h = z[0] - b;
        Array<2> res = Array<2>{ 0.5 * h * (z[1] * z[1] + z[2] * z[2]), 0.5 * h * h + h * b};
        return res;
      });
    }
  }

  return res;
}

Storage<3> Compare(const SpaceDisc& sd, const Test& test, double t) {
  const auto& m = sd.GetDomain().GetTopology();
  const size_t n = m.NumTriangles();

  Storage<3> errL2;
  errL2.resize(Eigen::NoChange, n);
  std::aligned_storage_t<sizeof(MUSCLObject::MUSCL), alignof(MUSCLObject::MUSCL)> data;
  MUSCLObject::MUSCL* muscl = nullptr;
  auto comparator = [&](const Point& p) {
    double b = sd.GetDomain().AtPoint(muscl->OriginId(), p);
    Array<3> U = Array<3>{ test.w(p[0], p[1], t), test.hu(p[0], p[1], t), test.hv(p[0], p[1], t) };
    Array<3> V = muscl->AtPoint(p);
    V.tail<2>() *= (V[0] - b);
    Array<3> res = (U - V).abs();
    return res;
  };

  for (size_t i = 0; i < n; ++i) {
    const auto& ip = m.TriangPoints(i);
    if (sd.IsFullWetCell(i)) {
      muscl = new (&data) MUSCLObject::MUSCL {sd.ReconstructFullWetCell(i)};
      errL2.col(i) = TriangAverage<3,1>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), comparator);
    } else if (sd.IsPartWetCell(i)) {
      muscl = new (&data) MUSCLObject::MUSCL {sd.ReconstructPartWetCell2(i)};
      errL2.col(i) = TriangAverage<3,100>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), comparator);
    } else {
      muscl = new (&data) MUSCLObject::MUSCL {sd.ReconstructDryCell(i)};
      errL2.col(i) = TriangAverage<3,1>(m.P(ip[0]), m.P(ip[1]), m.P(ip[2]), comparator);
    }
 }

  return errL2;
}
*/
