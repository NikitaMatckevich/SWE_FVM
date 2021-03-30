#pragma once
#include <TriangMesh.h>
#include <PointOperations.h>

bool is_wet(double h) noexcept;
Eigen::Array3d dry_state(double b) noexcept;

struct Bathymetry {
  Bathymetry(TriangMesh&& m);
  Bathymetry(Bathymetry&& other) = default;
  inline TriangMesh const& mesh() const { return m_; }
  inline Storage1d  const& buff() const { return b_; }
  inline size_t str_size() const { return b_.size(); }

  double& at_node (idx n);

  double  at_point(idx t, Point const& p) const;
	double  at_node (idx n) const;
  double  at_edge (idx n) const;
	double  at_cell (idx n) const;

 private:
  TriangMesh m_;
  Storage1d  b_;
 public:
  auto    at_nodes(NodeTagArray const& ns) const -> decltype(b_.operator()(ns));
};

struct BaseBathymetryWrapper {
  BaseBathymetryWrapper(Bathymetry const& b);
 protected:
  Bathymetry const& b_;
};
