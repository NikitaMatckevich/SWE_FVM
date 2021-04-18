#pragma once
#include <PointOperations.h>
#include <TriangMesh.h>

bool IsWet(double h) noexcept;

Array<3> DryState(double b) noexcept;

struct Bathymetry {
  
  explicit Bathymetry(TriangMesh&& m);

  inline const TriangMesh& Mesh() const { return m_m; }
  inline const Storage<1>& Buff() const { return m_b; }
  inline size_t Size() const { return m_b.size(); }

  double& AtNode (Idx n);
  double  AtNode (Idx n) const;
  double  AtPoint(Idx t, const Point& p) const;
	
  Array<2> Gradient(Idx n) const;

 private:
  TriangMesh  m_m;
  Storage<1>  m_b;

 public:
  auto AtNodes(NodeTagArray const& ns) const -> decltype(m_b.operator()(ns));
};

struct BaseBathymetryWrapper {
  explicit BaseBathymetryWrapper(Bathymetry const& b);
 protected:
  Bathymetry const& m_b;
};
