#include <TriangMesh.h>

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
        const Eigen::Ref<NodeTagArray>& nodes,
        const Eigen::Ref<EdgeTagArray>& edgeNodes,
        const Eigen::Ref<EdgeTagArray>& edgeElements,
        const Eigen::Ref<TriangTagArray>& elementNodes,
        const Eigen::Ref<TriangTagArray>& elementEdges,
        const Eigen::Ref<TriangTagArray>& elementNeighbours)
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

