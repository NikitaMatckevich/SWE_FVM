#include "Bathymetry.h"
#include "PointOperations.h"
using namespace Eigen;

Bathymetry::Bathymetry(TriangMesh const& M, ArrayXd const& z) : M_(M), z_(z) {}

double Bathymetry::operator()(index t, Array2d const& p) const {
  triang_tag const& tp = M_.triang_points(t);
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
double Bathymetry::p(index p) const {
  return z_[p];
}
double Bathymetry::t(index t) const {
  triang_tag const& tp = M_.triang_points(t);
  return 1. / 3. * (z_[tp[0]] + z_[tp[1]] + z_[tp[2]]);
}
double Bathymetry::e(index e) const {
  edge_tag const& ep = M_.edge_points(e);
  return 0.5 * (z_[ep[0]] + z_[ep[1]]);
}
double Bathymetry::c(index e) const {
  edge_tag const& ep = M_.edge_points(e);
  double l1 = length(M_.p(ep[0]), M_.c(e));
  double l2 = length(M_.p(ep[1]), M_.c(e));
  return (l2 * z_[ep[0]] + l1 * z_[ep[1]]) / (l1 + l2);
}
ArrayXd Bathymetry::p(node_array const& i) const {
  return z_(i);
}