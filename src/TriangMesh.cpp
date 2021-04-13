#include <TriangMesh.h>
#include <string>
#include <fstream>
#include <functional>
using namespace std;

// PARSING OF MSH FILE

TriangMesh::TriangMesh(string const& filename) {
  ifstream msh(filename);
  if (!msh.is_open())
    throw MeshError("cannot open mesh file " + filename);
  string line;
  auto find_section = bind([&filename](ifstream& msh, string& line, string const& sectname) {
    while (getline(msh, line) && line.find(sectname) != 0);
    if (msh.eof())
      throw MeshError("EOF while looking for " + sectname + " section in " + filename);
  }, std::ref(msh), std::ref(line), std::placeholders::_1);
  auto check_tag_format_too_small_to_handle = [&](idx i) {
    if (i >= numeric_limits<idx>::max())
      throw MeshError("idx is too small to handle in " + filename);
  };
  constexpr auto msh_size = numeric_limits<streamsize>::max();

  find_section("$MeshFormat");
  double version;
  idx is_binary;
  msh >> version >> is_binary;
  if (version != 4.0)
    throw MeshError("mesh format in " + filename + " is not of 4.0 version, the current version is: " + to_string(version));
  if (is_binary == 1)
    throw MeshError("mesh format in " + filename + " is binary, not supported");
  msh.ignore(msh_size, '\n');

  find_section("$Entities");
  idx num_entity_points, num_entity_curves;
  msh >> num_entity_points >> num_entity_curves;
  check_tag_format_too_small_to_handle(num_entity_curves);
  msh.ignore(msh_size, '\n');
  while (getline(msh, line) && (num_entity_points-- > 0));
  unordered_map<idx, idx> boundary_types;
  while (num_entity_curves-- > 0) {
    idx tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    idx num_phys_tags, boundary_type;
    msh >> tag >>
      min_x >> min_y >> min_z >>
      max_x >> max_y >> max_z >>
      num_phys_tags;
    if (num_phys_tags == 1) {
      msh >> boundary_type;
      boundary_types[tag] = -boundary_type;
    }
    msh.ignore(msh_size, '\n');
  }

  find_section("$Nodes");
  idx num_entity_blocks, num_nodes;
  msh >> num_entity_blocks >> num_nodes;
  check_tag_format_too_small_to_handle(num_entity_blocks);
  check_tag_format_too_small_to_handle(num_nodes);
  p_.resize(Eigen::NoChange, num_nodes);
  idx current_node = 1;
  for (idx block = 1; block <= num_entity_blocks; block++) {
    idx entity_tag, entity_dim, parametric, nodes_in_block;
    msh >> entity_tag >> entity_dim >> parametric >> nodes_in_block;
    for (idx node = 1; node <= nodes_in_block; node++) {
      idx ind;
      double x, y;
      msh >> ind >> x >> y;
      if (ind != current_node++)
        throw MeshError("node tags in " + filename + " are not consequent");
      p_(0, current_node) = x;
      p_(1, current_node) = y;
      msh.ignore(msh_size, '\n');
    }
  }

  find_section("$Elements");
  idx num_elements;
  msh >> num_entity_blocks >> num_elements;
  check_tag_format_too_small_to_handle(num_entity_blocks);
  check_tag_format_too_small_to_handle(num_elements);
}

// EIGEN PUBLIC INTERFACE

using namespace Eigen;
#include "PointOperations.h"

using NodeTag = TriangMesh::NodeTag;

EdgeTag TriangMesh::edge_points(NodeTag i) const { return edge_points_.row(i); }
EdgeTag TriangMesh::edge_triangs(NodeTag i) const { return edge_triangs_.row(i); }
bool TriangMesh::is_edge_boundary(NodeTag i) const {
  return edge_triangs_(i, 1) < 0;
}

TriangTag TriangMesh::triang_points(NodeTag i) const {
  return triang_points_.row(i);
}
TriangTag TriangMesh::triang_edges(NodeTag i) const {
  return triang_edges_.row(i);
}
TriangTag TriangMesh::triang_triangs(NodeTag i) const {
  return triang_triangs_.row(i);
}
bool TriangMesh::is_triangle_boundary(NodeTag i) const {
  auto ie = triang_edges_.row(i);
  return ie.unaryExpr([this](NodeTag i) {
    return this->is_edge_boundary(i);
  }).any();
}

Point TriangMesh::p(NodeTag i) const {
  return p_.col(i);
}
Point TriangMesh::t(NodeTag i) const {
  auto const& ip = triang_points_.row(i);
  return 1. / 3. * (p_.col(ip[0]) + p_.col(ip[1]) + p_.col(ip[2]));
}
Point TriangMesh::e(NodeTag i) const {
  auto const& ip = edge_points_.row(i);
  return 0.5 * (p_.col(ip[0]) + p_.col(ip[1]));
}
Point TriangMesh::c(NodeTag i) const {
  auto const& ip = edge_points_.row(i);
  if (is_edge_boundary(i)) return 0.5 * (p(ip[0]) + p(ip[1]));
  auto const& it = edge_triangs_.row(i);
  return intersection(p(ip[0]), p(ip[1]), t(it[0]), t(it[1]));
}
PointArray TriangMesh::p(NodeTagArray const& i) const {
  return p_(all, i);
}
PointArray TriangMesh::t(NodeTagArray const& i) const {
  PointArray res(2, i.rows());
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = t(i[k]);
  return res;
}
PointArray TriangMesh::e(NodeTagArray const& i) const {
  PointArray res(i.rows(), 1);
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = e(i[k]);
  return res;
}
PointArray TriangMesh::c(NodeTagArray const& i) const {
  PointArray res(i.rows(), 1);
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = c(i[k]);
  return res;
}

Vector2d TriangMesh::norm(NodeTag ie, NodeTag it) const {
  NodeTag a = edge_points_(ie, 0);
  NodeTag b = edge_points_(ie, 1);
  Vector2d tan = (p(b) - p(a)).matrix() / l(ie);
  if (det(t(it) - p(a), tan) < 0.) {
		tan = -tan;
	}
  return (Matrix2d() << 0, -1, 1, 0).finished() * tan;
}

double TriangMesh::l   (NodeTag ie) const {
  auto const& ip = edge_points_.row(ie);
  return len(p(ip[0]), p(ip[1]));
}

double TriangMesh::area(NodeTag it) const {
  auto const& ip = triang_points_.row(it);
  return triang_area(p(ip[0]), p(ip[1]), p(ip[2]));
}

double max_triang_area(TriangMesh const& M) {
  double max_S = 0.;
  for (size_t i = 0; i < M.num_triangles(); i++) {
    double S = M.area(i);
    if (S > max_S) max_S = S;
  }
  return max_S;
}

// SIMPLE I/O

#include <iostream>

void info(ostream& out, TriangMesh const& M) {
  out.setf(ios::left, ios::adjustfield);
  auto write_comment = bind(
    [&](ostream& out, idx wid, string_view c) {out.width(wid); out << c << ' '; },
    ref(out), std::placeholders::_1, std::placeholders::_2);
  auto write_line = [&](auto n, string_view c1, auto v1,
                                string_view c2, auto v2,
                                string_view c3, auto v3,
                                string_view c4, auto v4)
  {
    out.width(4); out << n  << ' ';
    write_comment(7, c1);
    out.width(3); out << v1 << ' ';
    write_comment(10, c2);
    out.width(3); out << v2 << ' ';
    write_comment(11, c3);
    out.width(3); out << v3 << ' ';
    write_comment(14, c4);
    out.width(3); out << v4 << ' ';
    out << '\n';
  };

  out << "TRIANGULAR MESH :\n";

  out << "nodes: total = " << M.num_nodes() << '\n';
  for (size_t i = 0; i < M.num_nodes(); ++i) {
    write_line(i, "points:", M.p(i).transpose(),
                  ""       , "",
                  ""       , "",
                  ""       , "");
  }
  
  out << "edges: total = " << M.num_edges() << '\n';
  for (size_t i = 0; i < M.num_edges(); ++i) {
    write_line(i, "points:"   , M.edge_points(i),
                  "triangles:", M.edge_triangs(i),
                  "center at:", M.e(i).transpose(),
                  ""          , "");
  }
  
  out << "triangles: total = " << M.num_triangles() 
      << ";\tmax area = " << max_triang_area(M) << '\n';
  for (size_t i = 0; i < M.num_triangles(); ++i) {
    write_line(i, "points:"       , M.triang_points(i),
                  "edges:"        , M.triang_edges(i),
                  "neighbours:"   , M.triang_triangs(i),
                  "barycenter at:", M.e(i).transpose());
  }
}
