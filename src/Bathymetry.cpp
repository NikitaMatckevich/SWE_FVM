#include <PointOperations.h>
#include <Bathymetry.h>

Point DryState(double b) noexcept { return { b, 0., 0. }; }

Domain::Domain(const Eigen::Ref<Storage<3>>& geometry, const Topology& topology)
        : m_geometry(geometry)
        , m_topology(topology) {}

Eigen::Vector2d Domain::TriangSlope(Idx t) const {
    const auto& tp = m_topology.TriangPoints(t);
    return Gradient(P(tp).matrix());
}

double Domain::Z(Idx t, double x, double y) const {
    const auto& tp = m_topology.TriangPoints(t);
	return TriangSlope(t).dot(Eigen::Vector2d(x, y) - P(tp[0]).topRows(2).matrix());
}

Point Domain::P(NodeTag i) const {
    return m_geometry.col(i);
}

Point Domain::T(NodeTag t) const {
    const auto& tp = m_topology.TriangPoints(t);
    return (P(tp).matrix() * Eigen::Vector3d::Constant(1. / 3.)).array();
}

Point Domain::E(NodeTag e) const {
    return Point{};
}

Point Domain::C(NodeTag e) const {
    if (m_topology.IsEdgeBoundary(e)) {
        return E(e);
    }
    const auto& ep = m_topology.EdgePoints(e);
    const auto& et = m_topology.EdgeTriangs(e);
    return Intersection(P(ep[0]), P(ep[1]), T(et[0]), T(et[1]));
}

PointArray Domain::T(const NodeTagArray& t) const {
    PointArray res;
    res.resize(Eigen::NoChange, t.size());
    for (NodeTag k = 0; k < t.size(); k++) {
        res.col(k) = T(t[k]);
    }
    return res;
}

PointArray Domain::E(const NodeTagArray& e) const {
    PointArray res;
    res.resize(Eigen::NoChange, e.size());
    for (NodeTag k = 0; k < e.size(); k++) {
        res.col(k) = E(e[k]);
    }
    return res;
}

PointArray Domain::C(const NodeTagArray& e) const {
  PointArray res;
  res.resize(Eigen::NoChange, e.size());
  for (NodeTag k = 0; k < e.size(); k++) {
    res.col(k) = C(e[k]);
  }
  return res;
}

Eigen::Vector2d Domain::Tang(NodeTag e, NodeTag t) const {
    const auto& ep = m_topology.EdgePoints(e);
    Eigen::Vector3d tan = (P(ep[1]) - P(ep[0])).matrix() / L(e);
    if (Det(T(t) - P(ep[0]), tan) > 0.) {
		tan = -tan;
    }
    return tan.topRows(2);
}

Eigen::Vector2d Domain::Norm(NodeTag e, NodeTag t) const {
  return (Eigen::Matrix2d() << 0, 1, -1, 0).finished() * Tang(e, t);
}

double Domain::L(NodeTag ie) const {
  const auto& ip = m_topology.EdgePoints(ie);
  return Len(P(ip[0]), P(ip[1]));
}

double Domain::Area(NodeTag it) const {
  const auto& ip = m_topology.TriangPoints(it);
  return TriangArea(P(ip[0]), P(ip[1]), P(ip[2]));
}

double MaxTriangArea(const Domain& b) {
    const auto& m = b.GetTopology();
    double maxArea = 0.;
    for (int i = 0; i < m.NumTriangles(); i++) {
        double area = b.Area(i);
        if (area > maxArea) {
            maxArea = area;
        }
    }
    return maxArea;
}

BaseDomainWrapper::BaseDomainWrapper(const Domain& b) : m_b(b) {}


// SIMPLE I/O

#include <iostream>
#include <fstream>
#include <functional>
#include <string>

using namespace std;

void Info(ostream& out, const Domain& b) {

    const auto& m = b.GetTopology();

    out.setf(ios::left, ios::adjustfield);

    auto WriteComment = bind(
        [&](ostream& out, Idx wid, string_view c) { out.width(wid); out << c << ' '; },
        ref(out), std::placeholders::_1, std::placeholders::_2);
    auto WriteLine = [&](auto n, string_view c1, auto v1,
                               string_view c2, auto v2,
                               string_view c3, auto v3,
                               string_view c4, auto v4)
    {
    out.width(4); out << n  << ' ';
    WriteComment(7,  c1);
    out.width(3); out << v1 << ' ';
    WriteComment(10, c2);
    out.width(3); out << v2 << ' ';
    WriteComment(11, c3);
    out.width(3); out << v3 << ' ';
    WriteComment(14, c4);
    out.width(3); out << v4 << ' ';
    out << '\n';
    };

    out << "TRIANGULAR MESH :\n";

    out << "nodes: total = " << m.NumNodes() << '\n';
    for (size_t i = 0; i < m.NumNodes(); ++i) {
        WriteLine(i, "points:", b.P(i).transpose(),
            "", "",
            "", "",
            "", "");
    }
  
    out << "edges: total = " << m.NumEdges() << '\n';
    for (size_t i = 0; i < m.NumEdges(); ++i) {
        WriteLine(i, "points:", m.EdgePoints(i),
            "triangles:", m.EdgeTriangs(i),
            "center at:", b.E(i).transpose(),
            "", "");
    }
  
    out << "triangles: total = " << m.NumTriangles() 
      << ";\tmax area = " << MaxTriangArea(b) << '\n';
    for (size_t i = 0; i < m.NumTriangles(); ++i) {
        WriteLine(i, "points:", m.TriangPoints(i),
            "edges:", m.TriangEdges(i),
            "nieghbours:", m.TriangTriangs(i),
            "center at:", b.T(i).transpose());
    }
}









