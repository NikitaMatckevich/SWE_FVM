#include <TriangMesh.h>
 
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

