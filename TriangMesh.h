#pragma once
#include "PointOperations.h"

using EdgeTag = Eigen::Array<index, 1, 2, Eigen::RowMajor>;
using TriangTag = Eigen::Array<index, 1, 3, Eigen::RowMajor>;

using NodeTagArray = Eigen::Array<index, 1, Eigen::Dynamic, Eigen::RowMajor>;
using EdgeTagArray = Eigen::Array<index, Eigen::Dynamic, 2, Eigen::RowMajor>;
using TriangTagArray = Eigen::Array<index, Eigen::Dynamic, 3, Eigen::RowMajor>;

struct MeshError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct TriangMesh {
  using NodeTag = index;
  TriangMesh() = default; //TODO: delete this
  TriangMesh(std::string const& filename);

  // TOPOLOGY SECTION
  inline size_t num_nodes() const { return p_.row(0).size(); }
  inline size_t num_edges() const { return edge_points_.col(0).size(); }
  inline size_t num_triangles() const { return triang_points_.col(0).size(); }

  EdgeTag edge_points(NodeTag i) const;
  EdgeTag edge_triangs(NodeTag i) const;
  bool is_edge_boundary(NodeTag i) const;

  TriangTag triang_points(NodeTag i) const;
  TriangTag triang_edges(NodeTag i) const;
  TriangTag triang_triangs(NodeTag i) const;
  bool is_triangle_boundary(NodeTag i) const;


  // GEOMETRY SECTION
  inline double min_x() const { return p_.row(0).minCoeff(); }
  inline double max_x() const { return p_.row(0).maxCoeff(); }
  inline double min_y() const { return p_.row(1).minCoeff(); }
  inline double max_y() const { return p_.row(1).maxCoeff(); }

  Point p(NodeTag i) const;
  Point t(NodeTag i) const;
  Point e(NodeTag i) const;
  Point c(NodeTag i) const;
  PointArray p(NodeTagArray const& i) const;
  PointArray t(NodeTagArray const& i) const;
  PointArray e(NodeTagArray const& i) const;
  PointArray c(NodeTagArray const& i) const;

  Eigen::Vector2d norm(NodeTag ie, NodeTag it) const;
  double l(NodeTag ie) const;
  double area(NodeTag it) const;
protected: //TODO : make private
  PointArray     p_;              // vertices of mesh elements

  EdgeTagArray   edge_points_;    // ends of each edge

  EdgeTagArray   edge_triangs_;   // triangles that have edge in common
                                  // some special types are used for boundary
                                  // edge's "ghost cells"

  TriangTagArray triang_points_;  // points of each triangle
  TriangTagArray triang_edges_;   // edges of each triangle
  TriangTagArray triang_triangs_; // triangles of each triangle
};

double max_triang_area(TriangMesh const&);
void info(std::ostream&, TriangMesh const&);