#pragma once
#include <Exceptions.h>
#include <PointOperations.h>

using EdgeTag   = Eigen::Array<Idx, 1, 2, Eigen::RowMajor>;
using TriangTag = Eigen::Array<Idx, 1, 3, Eigen::RowMajor>;

using NodeTagArray   = Eigen::Array<Idx, 1, Eigen::Dynamic, Eigen::RowMajor>;
using EdgeTagArray   = Eigen::Array<Idx, Eigen::Dynamic, 2, Eigen::RowMajor>;
using TriangTagArray = Eigen::Array<Idx, Eigen::Dynamic, 3, Eigen::RowMajor>;

struct TriangMesh {
  using NodeTag = Idx;
  TriangMesh() = default; // TODO(nikitamatckevich): delete this
  TriangMesh(TriangMesh&& other) = default;
  TriangMesh(const TriangMesh& other) = delete;
  TriangMesh(std::string const& filename);

  // TOPOLOGY SECTION
  inline size_t NumNodes() const { return m_p.row(0).size(); }
  inline size_t NumEdges() const { return m_edge_points.col(0).size(); }
  inline size_t NumTriangles() const { return m_triang_points.col(0).size(); }

  EdgeTag EdgePoints(NodeTag i) const;
  EdgeTag EdgeTriangs(NodeTag i) const;
  bool IsEdgeBoundary(NodeTag i) const;

  TriangTag TriangPoints(NodeTag i) const;
  TriangTag TriangEdges(NodeTag i) const;
  TriangTag TriangTriangs(NodeTag i) const;
  bool IsTriangleBoundary(NodeTag i) const;

  // GEOMETRY SECTION
  inline double MinX() const { return m_p.row(0).minCoeff(); }
  inline double MaxX() const { return m_p.row(0).maxCoeff(); }
  inline double MinY() const { return m_p.row(1).minCoeff(); }
  inline double MaxY() const { return m_p.row(1).maxCoeff(); }

  Point P(NodeTag i) const;
  Point T(NodeTag i) const;
  Point E(NodeTag i) const;
  Point C(NodeTag i) const;
  PointArray P(NodeTagArray const& i) const;
  PointArray T(NodeTagArray const& i) const;
  PointArray E(NodeTagArray const& i) const;
  PointArray C(NodeTagArray const& i) const;

  Eigen::Vector2d Norm(NodeTag ie, NodeTag it) const;
  double L(NodeTag ie) const;
  double Area(NodeTag it) const;
 
 protected: // TODO(nikitamatckevich): make private
  PointArray     m_p;              // vertices of mesh elements

  EdgeTagArray   m_edge_points;    // ends of each edge

  EdgeTagArray   m_edge_triangs;   // triangles that have edge in common
                                  // some special types are used for boundary
                                  // edge's "ghost cells"

  TriangTagArray m_triang_points;  // points of each triangle
  TriangTagArray m_triang_edges;   // edges of each triangle
  TriangTagArray m_triang_triangs; // triangles of each triangle
};

double MaxTriangArea(TriangMesh const&);

void Info(std::ostream&, TriangMesh const&);
