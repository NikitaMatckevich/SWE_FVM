#include <TriangMesh.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>

using namespace std;

struct EdgeHash {
  std::size_t operator()(const EdgeTag& e) const {
    return std::hash<Idx>()(e[0]) ^ (std::hash<Idx>()(e[1]) << 1);
 }
};
 
struct EdgeEqual {
  bool operator()(const EdgeTag& lhs, const EdgeTag& rhs) const {
    return lhs[0] == rhs[0] && lhs[1] == rhs[1];  
  }
};

Topology::Topology(
        const Eigen::Ref<const NodeTagArray>& nodes,
        const Eigen::Ref<const EdgeTagArray>& edgeNodes,
        const Eigen::Ref<const EdgeTagArray>& edgeElements,
        const Eigen::Ref<const TriangTagArray>& elementNodes,
        const Eigen::Ref<const TriangTagArray>& elementEdges,
        const Eigen::Ref<const TriangTagArray>& elementNeighbours)
    : m_nodes(nodes)
    , m_edge_points(edgeNodes)
    , m_edge_triangs(edgeElements)
    , m_triang_points(elementNodes)
    , m_triang_edges(elementEdges)
    , m_triang_triangs(elementNeighbours) {}

using namespace Eigen;

#include "PointOperations.h"

using NodeTag = Topology::NodeTag;

EdgeTag Topology::EdgePoints(NodeTag i) const { return m_edge_points.row(i); }
EdgeTag Topology::EdgeTriangs(NodeTag i) const { return m_edge_triangs.row(i); }
bool Topology::IsEdgeBoundary(NodeTag i) const {
  return m_edge_triangs(i, 1) < 0;
}

TriangTag Topology::TriangPoints(NodeTag i) const {
  return m_triang_points.row(i);
}
TriangTag Topology::TriangEdges(NodeTag i) const {
  return m_triang_edges.row(i);
}
TriangTag Topology::TriangTriangs(NodeTag i) const {
  return m_triang_triangs.row(i);
}
bool Topology::IsTriangleBoundary(NodeTag i) const {
  auto ie = m_triang_edges.row(i);
  return ie.unaryExpr([this](NodeTag i) {
    return this->IsEdgeBoundary(i);
  }).any();
}

double MaxTriangArea(Topology const& m) {
  double max_s = 0.;
  for (size_t i = 0; i < m.NumTriangles(); i++) {
    double s = m.Area(i);
    if (s > max_s) max_s = s;
  }
  return max_s;
}

// SIMPLE I/O

#include <iostream>

void Info(ostream& out, const Topology& m) {

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
    WriteLine(i, "points:", m.P(i).transpose(),
                  ""      , "",
                  ""      , "",
                  ""      , "");
  }
  
  out << "edges: total = " << m.NumEdges() << '\n';
  for (size_t i = 0; i < m.NumEdges(); ++i) {
    WriteLine(i, "points:"    , m.EdgePoints(i),
                  "triangles:", m.EdgeTriangs(i),
                  "center at:", m.E(i).transpose(),
                  ""          , "");
  }
  
  out << "triangles: total = " << m.NumTriangles() 
      << ";\tmax area = " << MaxTriangArea(m) << '\n';
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
    WriteLine(i, "points:"        , m.TriangPoints(i),
                  "edges:"        , m.TriangEdges(i),
                  "neighbours:"   , m.TriangTriangs(i),
                  "barycenter at:", m.E(i).transpose());
  }
}
