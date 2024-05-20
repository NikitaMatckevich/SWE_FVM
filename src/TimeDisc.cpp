#include <TimeDisc.h>

Array<3> TimeDisc::RHS(Idx i, double dt) const {

  const auto& vol = m_sd->GetVolField();
  const auto& edg = m_sd->GetEdgField();
  const auto& src = m_sd->GetSrcField();
  const auto& f   = m_sd->GetFluxes();

  Array<3> res = Array<3>::Zero();

  const auto& m  = m_sd->Mesh();
  const auto& ie = m.TriangEdges(i);
  const auto& it = m.TriangTriangs(i);
  
  double i_area = 1. / m.Area(i);

  double dti = ComputeDrainingDt(i); 

  for (int k = 0; k < 3; k++) {

    int sgn = m.EdgeTriangs(ie[k])[0] == i ? 1 : -1;
    
    double dtik = ComputeDrainingDt(it[k]);
    double dtk = (sgn * f(0, ie[k])) > 0. ? std::min(dt, dti) : std::min(dt, dtik);

    double c_ek = i_area * m.L(ie[k]);
    double h_ek = edg.h(ie[k], i, it[k]);
    
    res -= dtk * sgn * c_ek * f.col(ie[k]);

    // TODO(nikitamatckevich): if doesn't work, delete lines 51, uncomment line 54
    res.tail<2>() -= dt * (1./3.) * src.vel(ie[k], i, it[k]) * h_ek;
    res.tail<2>() += dtk * (m.Norm(ie[k], i) * c_ek * (0.5 * h_ek * h_ek)).array();
  }
  
  //res -= dt * s.prim(i) * vol.h(i);
  //res.tail<2>() -= dt * m_sd->Bath().Gradient(i) * vol.h(i);
  return res;
}

double TimeDisc::ComputeDrainingDt(Idx i) const {
  if (i < 0)
    return std::numeric_limits<double>::infinity(); 

  if (m_sd->IsDryCell(i))
    return 0.;

  const auto& m  = m_sd->Mesh();
  const auto& ie = m.TriangEdges(i);

  const auto& f = m_sd->GetFluxes();

  double sum = 0.;
  for (int k = 0; k < 3; k++) {
    Idx itk = m.EdgeTriangs(ie[k])[0];
    double f_ek = f(0, ie[k]);
    sum += std::max(0., (i == itk ? f_ek : -f_ek));
  }

  return sum > tol ? 
    m.Area(i) * m_sd->GetVolField().h(i) / sum :
    std::numeric_limits<double>::infinity();
}
