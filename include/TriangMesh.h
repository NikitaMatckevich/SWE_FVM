#pragma once
#include <Exceptions.h>
#include <PointOperations.h>

using EdgeTag   = Eigen::Array<Idx, 1, 2, Eigen::RowMajor>;
using TriangTag = Eigen::Array<Idx, 1, 3, Eigen::RowMajor>;

using NodeTagArray   = Eigen::Array<Idx, 1, Eigen::Dynamic, Eigen::RowMajor>;
using EdgeTagArray   = Eigen::Array<Idx, Eigen::Dynamic, 2, Eigen::RowMajor>;
using TriangTagArray = Eigen::Array<Idx, Eigen::Dynamic, 3, Eigen::RowMajor>;

struct Topology {
  using NodeTag = Idx;

  Topology(Topology&& other) = default;
  Topology(const Topology& other) = delete;

  Topology(
          const Eigen::Ref<NodeTagArray> nodes,
          const Eigen::Ref<EdgeTagArray> edgeNodes,
          const Eigen::Ref<EdgeTagArray> edgeElements,
          const Eigen::Ref<TriangTagArray> elementNodes,
          const Eigen::Ref<TriangTagArray> elementEdges,
          const Eigen::Ref<TriangTagArray> elementNeighbours);

  inline size_t NumNodes() const { return m_nodes.row(0).size(); }
  inline size_t NumEdges() const { return m_edge_points.col(0).size(); }
  inline size_t NumTriangles() const { return m_triang_points.col(0).size(); }

  EdgeTag EdgePoints(NodeTag i) const;
  EdgeTag EdgeTriangs(NodeTag i) const;
  bool IsEdgeBoundary(NodeTag i) const;

  TriangTag TriangPoints(NodeTag i) const;
  TriangTag TriangEdges(NodeTag i) const;
  TriangTag TriangTriangs(NodeTag i) const;
  bool IsTriangleBoundary(NodeTag i) const;

private:
  const NodeTagArray   m_nodes;         // mesh nodes
  const EdgeTagArray   m_edge_points;   // ends of each edge
  const EdgeTagArray   m_edge_triangs;  // triangle pairs that have 1 edge in common
                                        // some special types are used for boundary
                                        // edge's "ghost cells"

  const TriangTagArray m_triang_points;  // points of each triangle
  const TriangTagArray m_triang_edges;   // edges of each triangle
  const TriangTagArray m_triang_triangs; // triangles of each triangle
};

double MaxTriangArea(const Topology&);

void Info(std::ostream&, const Topology&);
