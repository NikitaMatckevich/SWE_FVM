#pragma once
#include <PointOperations.h>
#include <TriangMesh.h>

constexpr inline bool IsWet(double h) noexcept {
  constexpr double h_min = 1e-12; // minimal posible water depth of "wet" cell  
  return h > h_min;
}

Point DryState(double b) noexcept;

struct Domain {

    Domain(const Eigen::Ref<PointArray>& geometry, const Topology& topology);

    inline const Topology& GetTopology() const noexcept { return m_topology; }
    inline const Eigen::Ref<PointArray>& GetGeometry() const { return m_geometry; }
    inline size_t Size() const { return m_geometry.size(); }
   
    inline double MinX() const { return m_geometry.row(0).minCoeff(); }
    inline double MaxX() const { return m_geometry.row(0).maxCoeff(); }
    inline double MinY() const { return m_geometry.row(1).minCoeff(); }
    inline double MaxY() const { return m_geometry.row(1).maxCoeff(); }
    double Z(NodeTag t, double x, double y) const;
    
    Point P(NodeTag p) const;
    Point T(NodeTag t) const;
    Point E(NodeTag e) const;
    Point C(NodeTag c) const;
    PointArray T(const NodeTagArray& t) const;
    PointArray E(const NodeTagArray& e) const;
    PointArray C(const NodeTagArray& e) const;
  
    Eigen::Vector2d Tang(NodeTag e, NodeTag t) const;
    Eigen::Vector2d Norm(NodeTag e, NodeTag t) const;
    Eigen::Vector2d TriangSlope(NodeTag t) const;

    double L(NodeTag e) const;
    double Area(NodeTag t) const;
   
private:
    const Topology& m_topology;
    const Eigen::Ref<PointArray>& m_geometry;

public:
    auto P(const NodeTagArray& ns) const -> decltype(m_geometry.operator()(Eigen::all, ns)) {
        return m_geometry(Eigen::all, ns);
    }
};

struct BaseDomainWrapper {
    explicit BaseDomainWrapper(const Domain& b);
protected:
    const Domain& m_b;
};

double MaxTriangArea(const Domain&);

void Info(std::ostream&, const Domain&);
