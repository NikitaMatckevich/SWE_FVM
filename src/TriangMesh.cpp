#include <TriangMesh.h>
#include <fstream>
#include <functional>
#include <string>
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
 
  auto check_tag_format_too_small_to_handle = [&](Idx i) {
    if (i >= numeric_limits<Idx>::max())
      throw MeshError("Idx is too small to handle in " + filename);
  };

  constexpr auto msh_size = numeric_limits<streamsize>::max();

  find_section("$MeshFormat");
  double version;
  Idx is_binary;
  msh >> version >> is_binary;
  if (version != 4.0)
    throw MeshError("mesh format in " + filename + " is not of 4.0 version, the current version is: " + to_string(version));
  if (is_binary == 1)
    throw MeshError("mesh format in " + filename + " is binary, not supported");
  msh.ignore(msh_size, '\n');

  find_section("$Entities");
  Idx num_entity_points, num_entity_curves;
  msh >> num_entity_points >> num_entity_curves;
  check_tag_format_too_small_to_handle(num_entity_curves);
  msh.ignore(msh_size, '\n');
  while (getline(msh, line) && (num_entity_points-- > 0));
  unordered_map<Idx, Idx> boundary_types;
  while (num_entity_curves-- > 0) {
    Idx tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    Idx num_phys_tags, boundary_type;
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
  Idx num_entity_blocks, num_nodes;
  msh >> num_entity_blocks >> num_nodes;
  check_tag_format_too_small_to_handle(num_entity_blocks);
  check_tag_format_too_small_to_handle(num_nodes);
  m_p.resize(Eigen::NoChange, num_nodes);
  Idx current_node = 1;
  for (Idx block = 1; block <= num_entity_blocks; block++) {
    Idx entity_tag, entity_dim, parametric, nodes_in_block;
    msh >> entity_tag >> entity_dim >> parametric >> nodes_in_block;
    for (Idx node = 1; node <= nodes_in_block; node++) {
      Idx ind;
      double x, y;
      msh >> ind >> x >> y;
      if (ind != current_node++)
        throw MeshError("node tags in " + filename + " are not consequent");
      m_p(0, current_node) = x;
      m_p(1, current_node) = y;
      msh.ignore(msh_size, '\n');
    }
  }

  find_section("$Elements");
  Idx num_elements;
  msh >> num_entity_blocks >> num_elements;
  check_tag_format_too_small_to_handle(num_entity_blocks);
  check_tag_format_too_small_to_handle(num_elements);
}

// EIGEN PUBLIC INTERFACE

using namespace Eigen;
#include "PointOperations.h"

using NodeTag = TriangMesh::NodeTag;

EdgeTag TriangMesh::EdgePoints(NodeTag i) const { return m_edge_points.row(i); }
EdgeTag TriangMesh::EdgeTriangs(NodeTag i) const { return m_edge_triangs.row(i); }
bool TriangMesh::IsEdgeBoundary(NodeTag i) const {
  return m_edge_triangs(i, 1) < 0;
}

TriangTag TriangMesh::TriangPoints(NodeTag i) const {
  return m_triang_points.row(i);
}
TriangTag TriangMesh::TriangEdges(NodeTag i) const {
  return m_triang_edges.row(i);
}
TriangTag TriangMesh::TriangTriangs(NodeTag i) const {
  return m_triang_triangs.row(i);
}
bool TriangMesh::IsTriangleBoundary(NodeTag i) const {
  auto ie = m_triang_edges.row(i);
  return ie.unaryExpr([this](NodeTag i) {
    return this->IsEdgeBoundary(i);
  }).any();
}

Point TriangMesh::P(NodeTag i) const {
  return m_p.col(i);
}
Point TriangMesh::T(NodeTag i) const {
  const auto& ip = m_triang_points.row(i);
  return 1. / 3. * (m_p.col(ip[0]) + m_p.col(ip[1]) + m_p.col(ip[2]));
}
Point TriangMesh::E(NodeTag i) const {
  const auto& ip = m_edge_points.row(i);
  return 0.5 * (m_p.col(ip[0]) + m_p.col(ip[1]));
}

Point TriangMesh::C(NodeTag i) const {
  const auto& ip = m_edge_points.row(i);
  if (IsEdgeBoundary(i)) return 0.5 * (P(ip[0]) + P(ip[1]));
  const auto& it = m_edge_triangs.row(i);
  return Intersection(P(ip[0]), P(ip[1]), T(it[0]), T(it[1]));
}

PointArray TriangMesh::P(NodeTagArray const& i) const {
  return m_p(all, i);
}

PointArray TriangMesh::T(NodeTagArray const& i) const {
  PointArray res(2, i.rows());
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = T(i[k]);
  return res;
}

PointArray TriangMesh::E(NodeTagArray const& i) const {
  PointArray res(i.rows(), 1);
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = E(i[k]);
  return res;
}

PointArray TriangMesh::C(NodeTagArray const& i) const {
  PointArray res(i.rows(), 1);
  for (NodeTag k = 0; k < i.size(); k++)
    res.col(k) = C(i[k]);
  return res;
}

Vector2d TriangMesh::Norm(NodeTag ie, NodeTag it) const {
  NodeTag a = m_edge_points(ie, 0);
  NodeTag b = m_edge_points(ie, 1);
  Vector2d tan = (P(b) - P(a)).matrix() / L(ie);
  if (Det(T(it) - P(a), tan) < 0.) {
		tan = -tan;
	}
  return (Matrix2d() << 0, -1, 1, 0).finished() * tan;
}

double TriangMesh::L(NodeTag ie) const {
  const auto& ip = m_edge_points.row(ie);
  return Len(P(ip[0]), P(ip[1]));
}

double TriangMesh::Area(NodeTag it) const {
  const auto& ip = m_triang_points.row(it);
  return TriangArea(P(ip[0]), P(ip[1]), P(ip[2]));
}

double MaxTriangArea(TriangMesh const& m) {

  double max_s = 0.;

  for (size_t i = 0; i < m.NumTriangles(); i++) {
    double s = m.Area(i);
    if (s > max_s) max_s = s;
  }
  return max_s;
}

// SIMPLE I/O

#include <iostream>

void Info(ostream& out, const TriangMesh& m) {

  out.setf(ios::left, ios::adjustfield);

  auto write_comment = bind(
    [&](ostream& out, Idx wid, string_view c) {out.width(wid); out << c << ' '; },
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

  out << "nodes: total = " << m.NumNodes() << '\n';
  for (size_t i = 0; i < m.NumNodes(); ++i) {
    write_line(i, "points:", m.P(i).transpose(),
                  ""       , "",
                  ""       , "",
                  ""       , "");
  }
  
  out << "edges: total = " << m.NumEdges() << '\n';
  for (size_t i = 0; i < m.NumEdges(); ++i) {
    write_line(i, "points:"   , m.EdgePoints(i),
                  "triangles:", m.EdgeTriangs(i),
                  "center at:", m.E(i).transpose(),
                  ""          , "");
  }
  
  out << "triangles: total = " << m.NumTriangles() 
      << ";\tmax area = " << MaxTriangArea(m) << '\n';
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
    write_line(i, "points:"       , m.TriangPoints(i),
                  "edges:"        , m.TriangEdges(i),
                  "neighbours:"   , m.TriangTriangs(i),
                  "barycenter at:", m.E(i).transpose());
  }
}
