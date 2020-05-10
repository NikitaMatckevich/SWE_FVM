#pragma once
#include "TriangMesh.h"
#include "PointOperations.h"

bool is_wet(double h) noexcept;
Eigen::Array3d dry_state(double b) noexcept;

struct Bathymetry {
  using ValueArray = Eigen::ArrayXd;
  Bathymetry(TriangMesh const& M, ValueArray const& z);
  inline TriangMesh const& mesh() const { return M_; }
  inline ValueArray const& buffer() const { return z_; }
  inline size_t size() const { return z_.size(); }

  double operator()(index t, Point const& p) const;
  double p(index p) const;
  double t(index t) const;
  double e(index e) const;
  double c(index e) const;
  ValueArray p(NodeTagArray const& i) const;
private:
  TriangMesh const& M_;
  ValueArray const& z_;
};