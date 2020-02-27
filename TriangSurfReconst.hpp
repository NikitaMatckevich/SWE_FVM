#pragma once
#include "TriangMesh.hpp"

class TriangSurfReconst {
private:
  TriangMesh const & M_;
  Eigen::ArrayXd const & z_;
public:
  TriangSurfReconst(TriangMesh const & M, Eigen::ArrayXd const & z) : M_(M), z_(z) {}
  TriangMesh const & mesh() const { return M_; }
  Eigen::ArrayXd const & buffer() const { return z_; }
  size_t size() const { return z_.size(); }

  double operator()(node_tag t, point const & p) const {
    triang_tag const & tp = M_.triang_points(t);
    double idet = 1. / det(M_.p(tp[2]) - M_.p(tp[0]), M_.p(tp[2]) - M_.p(tp[1]));
    double lam0 = det(M_.p(tp[2]) - p, M_.p(tp[2]) - M_.p(tp[1])) * idet;
    assert(lam0 >= 0.);
    double lam1 = det(p - M_.p(tp[2]), M_.p(tp[2]) - M_.p(tp[0])) * idet;
    assert(lam1 >= 0.);
    double lam2 = 1. - lam0 - lam1;
    assert(lam2 >= 0.);
    return z_[tp[0]] * lam0 +
      z_[tp[1]] * lam1 +
      z_[tp[2]] * lam2;
  }
  double p(node_tag p) const {
    return z_[p];
  }
  double t(node_tag t) const {
    triang_tag const & tp = M_.triang_points(t);
    return 1. / 3. * (z_[tp[0]] + z_[tp[1]] + z_[tp[2]]);
  }
  double e(node_tag e) const {
    edge_tag const & ep = M_.edge_points(e);
    return 0.5 * (z_[ep[0]] + z_[ep[1]]);
  }
  double c(node_tag e) const {
    edge_tag const & ep = M_.edge_points(e);
    double l1 = length(M_.p(ep[0]), M_.c(e));
    double l2 = length(M_.p(ep[1]), M_.c(e));
    return (l2 * z_[ep[0]] + l1 * z_[ep[1]]) / (l1 + l2);
  }
  ArrayXd p(node_array const & i) const {
    return z_(i);
  }
};

double avg_triang(point const & p0, point const & p1, point const & p2,
  std::function<double(point const &)> const & f) {
  double S = triangle_area(p0, p1, p2);
  const unsigned int N = 500;
  double h = 1. / N;
  point di = h * (p1 - p0);
  point dj = h * (p2 - p0);
  point dt = 1. / 3. * (di + dj);
  double sum = 0.;
  point pi = p0;
  for (unsigned int i = 0; i < N; i++) {
    point pt = pi + dt;
    for (unsigned int j = 0; j < N - i - 1; j++) {
      sum += h * f(pt);
      sum += h * f(pt + dt);
      pt += dj;
    }
    sum += h * f(pt);
    pi += di;
  }
  return h * sum;
}