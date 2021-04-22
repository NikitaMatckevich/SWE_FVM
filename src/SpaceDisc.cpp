#include <SpaceDisc.h>

void SpaceDisc::ComputeFluxes() {

  m_min_length_to_wavespeed = 1.;
  const auto& m = Mesh();

  for (size_t i = 0; i < m.NumEdges(); ++i) {

    const auto& it = m.EdgeTriangs(i);
    Idx lf = it[0];
    Idx lt = it[1];

    //if (IsDryCell(lf) && IsDryCell(lt)) {
    //  m_f.col(i) = Array<3>{0., 0., 0.};
    //  continue;
    //}

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

void SpaceDisc::DumpFluxes(const std::string& filename) const {

  std::ofstream fout(filename);

  fout << "x\ty\tf1\tf2\tf3\tt1\tt2\n";

  for (size_t i = 0; i < Mesh().NumEdges(); i++) {
    
    const auto& e = Mesh().E(i);
    const auto& it = Mesh().EdgeTriangs(i);
    
    fout << e[0] << '\t' << e[1] << '\t';
    fout << m_f(0, i) << '\t' << m_f(1, i) << '\t' << m_f(2, i) << '\t';
    fout << it[0] << '\t' << it[1] << '\n';
  }

  fout.close();
}
