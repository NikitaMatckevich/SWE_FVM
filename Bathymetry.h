#pragma once
#include "TriangMesh.h"

class Bathymetry {
  TriangMesh const& M_;
  Eigen::ArrayXd const& z_;
public:
  Bathymetry(TriangMesh const& M, Eigen::ArrayXd const& z);
  inline TriangMesh const& mesh() const { return M_; }
  inline Eigen::ArrayXd const& buffer() const { return z_; }
  inline size_t size() const { return z_.size(); }

  double operator()(index t, Eigen::Array2d const& p) const;
  double p(index p) const;
  double t(index t) const;
  double e(index e) const;
  double c(index e) const;
  Eigen::ArrayXd p(node_array const& i) const;
};