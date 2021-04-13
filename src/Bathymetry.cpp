#include <Bathymetry.h>
#include <PointOperations.h>

bool is_wet(double h) noexcept {
  constexpr double h_min = 1e-10; // minimal posible water depth of "wet" cell  
  return h > h_min;
}

Array<3> dry_state(double b) noexcept { return { b, 0., 0. }; }

Bathymetry::Bathymetry(TriangMesh&& m) : m_(std::move(m)) {
  b_.resize(Eigen::NoChange, m_.num_nodes());
}
double  Bathymetry::at_point(idx t, Point const& p) const {
  auto const& tp = m_.triang_points(t);
  double idet = 1. / det(m_.p(tp[2]) - m_.p(tp[0]), m_.p(tp[2]) - m_.p(tp[1]));
  double lam0 = det(m_.p(tp[2]) - p, m_.p(tp[2]) - m_.p(tp[1])) * idet;
  assert(lam0 >= 0.);
  double lam1 = det(p - m_.p(tp[2]), m_.p(tp[2]) - m_.p(tp[0])) * idet;
  assert(lam1 >= 0.);
  double lam2 = 1. - lam0 - lam1;
  assert(lam2 >= 0.);
  return b_[tp[0]] * lam0 + b_[tp[1]] * lam1 + b_[tp[2]] * lam2;
}
double& Bathymetry::at_node(idx n) {
  return b_[n];
}
double  Bathymetry::at_node(idx n) const {
  return b_[n];
}
double  Bathymetry::at_edge(idx n) const {
	const auto& ep = mesh().edge_points(n); 
  return 0.5 * (b_[ep[0]] + b_[ep[1]]);
}
double  Bathymetry::at_cell(idx n) const {
  const auto& tp = mesh().triang_points(n);
	return (1./3.) * (b_[tp[0]] + b_[tp[1]] + b_[tp[2]]);
}

Array<2> Bathymetry::grad(idx n) const {
  const auto& tp = mesh().triang_points(n);
	return (at_nodes(tp).matrix() * gradient_coefs(mesh().p(tp).transpose())).array().transpose();
}

auto Bathymetry::at_nodes(NodeTagArray const& ns) const
-> decltype(b_.operator()(ns)) {
  return b_(ns);
}

BaseBathymetryWrapper::BaseBathymetryWrapper(Bathymetry const& b) : b_(b) {}
