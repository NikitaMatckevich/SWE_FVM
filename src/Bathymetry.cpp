#include <Bathymetry.h>
#include <PointOperations.h>

Array<3> DryState(double b) noexcept { return { b, 0., 0. }; }

Bathymetry::Bathymetry(const TriangMesh* m) : m_m(m) {
  m_b.resize(Eigen::NoChange, m_m->NumNodes());
}

double  Bathymetry::AtPoint(Idx t, Point const& p) const {
  auto const& tp = m_m->TriangPoints(t);
  double idet = 1. / Det(m_m->P(tp[2]) - m_m->P(tp[0]), m_m->P(tp[2]) - m_m->P(tp[1]));
  double lam0 = Det(m_m->P(tp[2]) - p, m_m->P(tp[2]) - m_m->P(tp[1])) * idet;
  assert(lam0 >= 0.);
  double lam1 = Det(p - m_m->P(tp[2]), m_m->P(tp[2]) - m_m->P(tp[0])) * idet;
  assert(lam1 >= 0.);
  double lam2 = 1. - lam0 - lam1;
  assert(lam2 >= 0.);
  return m_b[tp[0]] * lam0 + m_b[tp[1]] * lam1 + m_b[tp[2]] * lam2;
}

double& Bathymetry::AtNode(Idx n) {
  return m_b[n];
}

double  Bathymetry::AtNode(Idx n) const {
  return m_b[n];
}

Array<2> Bathymetry::Gradient(Idx n) const {
  const auto& tp = m_m->TriangPoints(n);
	return (AtNodes(tp).matrix() * GradientCoefs(m_m->P(tp).transpose())).array().transpose();
}

auto Bathymetry::AtNodes(NodeTagArray const& ns) const -> decltype(m_b.operator()(ns))
{
  return m_b(ns);
}

BaseBathymetryWrapper::BaseBathymetryWrapper(const Bathymetry& b) : m_b(b) {}
