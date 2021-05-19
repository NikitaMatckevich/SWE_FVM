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

// PARSING OF MSH FILE

TriangMesh::TriangMesh(const string& filename) {

  ifstream msh(filename);
  if (!msh.is_open()) {
    throw MeshError("cannot open mesh file " + filename);
  }

  string line;

  auto FindSection = bind([&filename](ifstream& msh, string& line, const string& sectname) {
    while (getline(msh, line) && line.find(sectname) != 0);
    if (msh.eof())
      throw MeshError("EOF while looking for " + sectname + " section in " + filename);
  }, std::ref(msh), std::ref(line), std::placeholders::_1); 
  constexpr auto msh_size = numeric_limits<streamsize>::max();

  FindSection("$MeshFormat");
  double version;
  Idx is_binary;
  msh >> version >> is_binary;
  if (version < 4.1) {
    throw MeshError("mesh format in "
      + filename
      + " is not >= 4.1 version, the current version is: "
      + to_string(version));
  }
  if (is_binary == 1) {
    throw MeshError("mesh format in "
      + filename
      + " is binary, not supported");
  }

  FindSection("$Entities");
  Idx num_entity_points, num_entity_curves, num_entity_blocks;
  msh >> num_entity_points >> num_entity_curves;
  msh.ignore(msh_size, '\n');
  vector<Idx> physical_groups(num_entity_curves);

  while (getline(msh, line) && (--num_entity_points > 0));
  for (Idx curve = 1; curve <= num_entity_curves; curve++) {
    Idx curve_tag, num_physical_tags, physical_tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    msh >> curve_tag;
    if (curve != curve_tag) {
      throw MeshError("curve tags in " + filename + " are not consequent");
    }
    msh >> min_x >> min_y >> min_z;
    msh >> max_x >> max_y >> max_z;
    msh >> num_physical_tags;
    if (num_physical_tags != 1) {
      throw MeshError("ambiguous definition of boundary condition at file " + filename);
    }
    msh >> physical_tag;
    physical_groups[curve-1] = physical_tag;
    msh.ignore(msh_size, '\n');
  }

  FindSection("$Nodes");
  Idx num_nodes;
  msh >> num_entity_blocks >> num_nodes;
  msh.ignore(msh_size, '\n');
  m_p.resize(Eigen::NoChange, num_nodes);

  for (Idx block = 1, current_node = 0; block <= num_entity_blocks; block++) {
    Idx entity_dim, entity_tag, parametric, nodes_in_block;
    msh >> entity_dim >> entity_tag >> parametric >> nodes_in_block;
    Idx cnt = nodes_in_block;
    while (getline(msh, line) && (cnt-- > 0));
    for (Idx node = 1; node <= nodes_in_block; node++) {
      double x, y;
      msh >> x >> y; 
      m_p(0, current_node) = x;
      m_p(1, current_node) = y;
      current_node++;
      msh.ignore(msh_size, '\n');
    }
  }

  FindSection("$Elements");
  Idx num_elements;
  msh >> num_entity_blocks >> num_elements;
  msh.ignore(msh_size, '\n');
  
  unordered_map<EdgeTag, Idx, EdgeHash, EdgeEqual> edge_ids;
  m_edge_triangs.resize(2 * num_nodes - 1, Eigen::NoChange);

  for (Idx block = 1; block < num_entity_blocks; block++) {
    Idx entity_dim, entity_tag, elem_type, elems_in_block;
    msh >> entity_dim >> entity_tag >> elem_type >> elems_in_block;
    if (elem_type != 1) {
      throw MeshError("unknown elements before triangles in "
        + filename
        + ", expected only boundary edges");
    }
    for (Idx elem = 1; elem <= elems_in_block; elem++) {
      Idx id, beg, end;
      msh >> id >> beg >> end;
      if (--beg > --end) std::swap(beg, end);
      edge_ids.emplace(make_pair(EdgeTag{beg, end}, edge_ids.size()));
      m_edge_triangs(edge_ids.size() - 1, 1) = -physical_groups[entity_tag-1];
      msh.ignore(msh_size, '\n');
    }
  }

  Idx entity_dim, entity_tag, elem_type, num_triangs;
  msh >> entity_dim >> entity_tag >> elem_type >> num_triangs;
  if (elem_type != 2) {
    throw MeshError("last element block of "
      + filename
      + " should be a block of triangles");
  }

  m_edge_triangs.conservativeResize(num_nodes + num_triangs - 1, Eigen::NoChange);
  m_triang_points.resize(num_triangs, Eigen::NoChange);
  m_triang_edges.resize(num_triangs, Eigen::NoChange);
  m_triang_triangs.resize(num_triangs, Eigen::NoChange);

  for (Idx it = 0; it < num_triangs; it++) {
    Idx elem_id, pts[3];
    msh >> elem_id >> pts[0] >> pts[1] >> pts[2];
    for (short int k = 0; k < 3; k++) {
      Idx pa = pts[k] - 1;
      Idx pb = pts[(k + 1) % 3] - 1;
      m_triang_points(it, k) = pa;
      if (pa > pb) std::swap(pa, pb);
      auto insertion = edge_ids.emplace(EdgeTag{pa, pb}, edge_ids.size());
      Idx edge_id = (*insertion.first).second;
      m_triang_edges(it, k) = edge_id;
      if (insertion.second) {
        m_edge_triangs(edge_id, 1) = it; 
      } else {
        m_edge_triangs(edge_id, 0) = it;
        Idx itk = m_edge_triangs(edge_id, 1);
        m_triang_triangs(it, k) = itk;
        if (itk >= 0) {
          Idx l;
          (m_triang_edges.row(itk) == edge_id).maxCoeff(&l);
          m_triang_triangs(itk, l) = it;
        }
      }
    }
  }

  m_edge_points.resize(edge_ids.size(), Eigen::NoChange);
  for (const auto& edge : edge_ids) {
    m_edge_points(edge.second, 0) = edge.first[0];
    m_edge_points(edge.second, 1) = edge.first[1];
  }
}

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

Vector2d TriangMesh::Tang(NodeTag ie, NodeTag it) const {
  NodeTag a = m_edge_points(ie, 0);
  NodeTag b = m_edge_points(ie, 1);
  Vector2d tan = (P(b) - P(a)).matrix() / L(ie);
  if (Det(T(it) - P(a), tan) > 0.) {
		tan = -tan;
	}
  return tan;
}

Vector2d TriangMesh::Norm(NodeTag ie, NodeTag it) const {
  return (Matrix2d() << 0, 1, -1, 0).finished() * Tang(ie, it);
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
