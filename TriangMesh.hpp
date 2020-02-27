#pragma once
#include "Mesh.hpp"
#include "PointOperations.hpp"
using namespace Eigen;
using node_tag = int;
using edge_tag = Array<node_tag, 1, 2, RowMajor>;
using triang_tag = Array<node_tag, 1, 3, RowMajor>;
using node_array = Array<node_tag, 1, Dynamic, RowMajor>;

//Triangular mesh
class TriangMesh : public Mesh {
public:
  using edge_array = Array<node_tag, Dynamic, 2, RowMajor>;
  using triang_array = Array<node_tag, Dynamic, 3, RowMajor>;

  inline size_t num_nodes() const { return p_.row(0).size(); }
  inline size_t num_edges() const { return edge_points_.col(0).size(); }
  inline size_t num_triangles() const { return triang_points_.col(0).size(); }

  edge_tag edge_points(node_tag i) const { return edge_points_.row(i); }
  edge_tag edge_triangs(node_tag i) const { return edge_triangs_.row(i); }
  bool is_edge_boundary(node_tag i) const {
    return edge_triangs_(i, 1) < 0;
  }

  triang_tag triang_points(node_tag i) const {
    return triang_points_.row(i);
  }
  triang_tag triang_edges(node_tag i) const {
    return triang_edges_.row(i);
  }
  triang_tag triang_triangs(node_tag i) const {
    return triang_triangs_.row(i);
  }
  bool is_triangle_boundary(node_tag i) const {
    triang_tag const & e = triang_edges_.row(i);
    return e.unaryExpr([this](node_tag i) {
      return this->is_edge_boundary(i);
    }).any();
  }

  point p(node_tag i) const {
    return p_.col(i);
  }
  point t(node_tag i) const {
    triang_tag const & ip = triang_points_.row(i);
    return 1. / 3. * (p_.col(ip[0]) + p_.col(ip[1]) + p_.col(ip[2]));
  }
  point e(node_tag i) const {
    edge_tag const & ip = edge_points_.row(i);
    return 0.5 * (p_.col(ip[0]) + p_.col(ip[1]));
  }
  point c(node_tag i) const {
    edge_tag const & ip = edge_points_.row(i);
    if (is_edge_boundary(i)) return 0.5 * (p(ip[0]) + p(ip[1]));
    edge_tag const & it = edge_triangs_.row(i);
    return intersection(p(ip[0]), p(ip[1]), t(it[0]), t(it[1]));
  }
  Array2Xd p(node_array const & i) const {
    return p_(all, i);
  }
  Array2Xd t(node_array const & i) const {
    Array2Xd res(2, i.rows());
    for (node_tag k = 0; k < i.size(); k++)
      res.col(k) = t(i[k]);
    return res;
  }
  Array2Xd e(node_array const & i) const {
    Array2Xd res(i.rows(), 1);
    for (node_tag k = 0; k < i.size(); k++)
      res.col(k) = e(i[k]);
    return res;
  }
  Array2Xd c(node_array const & i) const {
    Array2Xd res(i.rows(), 1);
    for (node_tag k = 0; k < i.size(); k++)
      res.col(k) = c(i[k]);
    return res;
  }

  Vector2d normal(node_tag e) const {
    node_tag i = edge_triangs_(e, 0);
    node_tag a = edge_points_(e, 0);
    node_tag b = edge_points_(e, 1);
    point tan = (p(b) - p(a)) / l(e);
    Vector2d nor;
    if (det(t(i) - p(a), tan) < 0.) {
      nor[0] = tan[1];  nor[1] = -tan[0];
    }
    else {
      nor[0] = -tan[1];  nor[1] = tan[0];
    }
    return nor;
  }
  double l(node_tag e) const {
    edge_tag const & k = edge_points_.row(e);
    return length(p(k[0]), p(k[1]));
  }
  double area(node_tag t) const {
    triang_tag const & k = triang_points_.row(t);
    return triangle_area(p(k[0]), p(k[1]), p(k[2]));
  }
  double min_x() const override { return p_.row(0).minCoeff(); }
  double max_x() const override { return p_.row(0).maxCoeff(); }
  double min_y() const override { return p_.row(1).minCoeff(); }
  double max_y() const override { return p_.row(1).maxCoeff(); }

protected:
  TriangMesh() {}

  //Vertices of mesh elements
  Array2Xd p_;

  //Edge tags
  edge_array edge_points_;
  edge_array edge_triangs_;
  //Triangle tags
  triang_array triang_points_;
  triang_array triang_edges_;
  triang_array triang_triangs_;
};

double max_triang_area(TriangMesh const & M) {
  double max_S = 0.0;
  for (node_tag i = 0; i < M.num_triangles(); i++) {
    double S = M.area(i);
    if (S >= max_S) max_S = S;
  }
  return max_S;
}