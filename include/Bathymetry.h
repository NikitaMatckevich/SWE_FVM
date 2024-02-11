#pragma once
#include <PointOperations.h>
#include <TriangMesh.h>

consteval bool IsWet(double h) noexcept;

Array<3> DryState(double b) noexcept;

struct Bathymetry {
  
    explicit Bathymetry(TriangMesh* m);

    inline const TriangMesh& Mesh() const { return *m_m; }
    inline const Storage<1>& Buff() const { return m_b; }
    inline size_t Size() const { return m_b.size(); }

    double& AtNode (Idx n);
    double  AtNode (Idx n) const;
    double  AtPoint(Idx t, const Point& p) const;
	
    Array<2> Gradient(Idx n) const;

private:
    TriangMesh *const m_m;
    Storage<1> m_b;

public:
    auto AtNodes(const NodeTagArray& ns) const -> decltype(m_b.operator()(ns));
};

struct BaseBathymetryWrapper {
    explicit BaseBathymetryWrapper(const Bathymetry& b);
protected:
    const Bathymetry& m_b;
};
