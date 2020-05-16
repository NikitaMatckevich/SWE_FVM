#pragma once
#include "TriangMesh.h"
#include "PointOperations.h"

bool is_wet(double h) noexcept;
Eigen::Array3d dry_state(double b) noexcept;

struct Bathymetry {
  using Storage = Eigen::ArrayXd;
  Bathymetry(TriangMesh&& m);
  Bathymetry(Bathymetry&& other) = default;
  inline TriangMesh const& mesh() const { return m_; }
  inline Storage    const& buff() const { return b_; }
  inline size_t str_size() const { return b_.size(); }

  double  at_point(index t, Point const& p) const;
  double& at_node (index n);
  double  at_node (index n) const;
private:
  TriangMesh m_;
  Storage    b_;
public:
  auto    at_nodes(NodeTagArray const& ns) const -> decltype(b_.operator()(ns));
};

struct BaseBathymetryWrapper {
  BaseBathymetryWrapper(Bathymetry const& b);
protected:
  Bathymetry const& b_;
};