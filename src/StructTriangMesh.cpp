#include <StructTriangMesh.h>

StructTriangMesh::StructTriangMesh(size_t ni, size_t nj, double h)
	: m_ni(ni)
  , m_nj(nj)
	, m_num_verts((ni + 1)* (nj + 1))
	, m_num_sqrs(ni * nj)
{
	size_t num_triangs = 4 * m_num_sqrs;
	m_triang_points.resize(num_triangs, Eigen::NoChange);
	m_triang_edges.resize(num_triangs, Eigen::NoChange);
	m_triang_triangs.resize(num_triangs, Eigen::NoChange);

	size_t num_pts = m_num_verts + m_num_sqrs;
	m_p.resize(Eigen::NoChange, num_pts);
	
  size_t num_edges = num_triangs + num_pts - 1;
	m_edge_points.resize(num_edges, Eigen::NoChange);
	m_edge_triangs.resize(num_edges, Eigen::NoChange);
	
  // fill p-array with square's vertices
	for (Idx j = 0; j < m_nj + 1; j++) {
		for (Idx i = 0; i < m_ni + 1; i++) {
			m_p(0, (m_ni + 1) * j + i) = i * h;
			m_p(1, (m_ni + 1) * j + i) = j * h;
		}
	}

	// fill p-array with square's centers
	for (Idx j = 0; j < m_nj; j++) {
		for (Idx i = 0; i < m_ni; i++) {
			m_p(0, m_num_verts + m_ni * j + i) = (i + 0.5) * h;
			m_p(1, m_num_verts + m_ni * j + i) = (j + 0.5) * h;
		}
	}

	// determine triangles and edges
	for (Idx j = 0; j < m_ni; j++) {
		for (Idx i = 0; i < m_ni; i++) {
			TriangB(i, j); //-bottom triangle of this square
			TriangR(i, j); //-right  triangle of this square
			TriangT(i, j); //-top    triangle of this square
			TriangL(i, j); //-left   triangle of this square
		}
	}
}

void StructTriangMesh::TriangB(Idx i, Idx j) {
	Idx it = (j * m_ni + i) * 4 + 0;
	Idx it0 = (j == 0) ? static_cast<Idx>(Boundaries::SOLID_WALL) :
		4 * ((j - 1) * m_ni + i) + 2;
	Idx it1 = (j * m_ni + i) * 4 + 1;
	Idx it2 = (j * m_ni + i) * 4 + 3;
	m_triang_triangs.row(it) << it0, it1, it2;
	Idx ip0 = (m_ni + 1) * j + i;
	Idx ip1 = (m_ni + 1) * j + i + 1;
	Idx ip2 = m_num_verts + m_ni * j + i;
	m_triang_points.row(it) << ip0, ip1, ip2;
	if (j == 0) {
		m_edge_points.row(m_curr_e) << ip0, ip1;
		m_edge_triangs.row(m_curr_e++) << it, it0;
	}
	m_edge_points.row(m_curr_e) << ip1, ip2;
	m_edge_triangs.row(m_curr_e++) << it, it1;
	m_edge_points.row(m_curr_e) << ip0, ip2;
	m_edge_triangs.row(m_curr_e++) << it, it2;
	if (j == 0) {
		// we've already put all 3 edges in mesh
		m_triang_edges.row(it) << m_curr_e - 3, m_curr_e - 2, m_curr_e - 1;
		return;
	}
	if (j == 1) {
		// take bottom edge from previous step
		// (it is actually top edge of some top
		// triangle in row j - 1)
		m_triang_edges.row(it) << m_curr_e - m_ni * 7 + i + 2, m_curr_e - 2, m_curr_e - 1;
		return;
	}
	// in general we do same thing as
	// in case (j == 1), but as we
	// have put more bottom edges in mesh
	// then in case (j == 1), we should
	// now decrement m_curr_e in other way
	m_triang_edges.row(it) << m_curr_e - m_ni * 6 + 1, m_curr_e - 2, m_curr_e - 1;
}

void StructTriangMesh::TriangR(Idx i, Idx j) {
	Idx it = (j * m_ni + i) * 4 + 1;
	Idx it0 = (i == m_ni - 1) ? static_cast<Idx>(Boundaries::SOLID_WALL) :
		4 * (j * m_ni + i + 1) + 3;
	Idx it1 = (j * m_ni + i) * 4 + 2;
	Idx it2 = (j * m_ni + i) * 4 + 0;
	m_triang_triangs.row(it) << it0, it1, it2;
	Idx ip0 = (m_ni + 1) * j + i + 1;
	Idx ip1 = (m_ni + 1) * (j + 1) + i + 1;
	Idx ip2 = m_num_verts + m_ni * j  + i;
	m_triang_points.row(it) << ip0, ip1, ip2;
	m_edge_points.row(m_curr_e) << ip0, ip1;
	m_edge_triangs.row(m_curr_e++) << it, it0;
	m_edge_points.row(m_curr_e) << ip1, ip2;
	m_edge_triangs.row(m_curr_e++) << it, it1;
	m_triang_edges.row(it) << m_curr_e - 2, m_curr_e - 1, m_curr_e - 4;
}

void StructTriangMesh::TriangT(Idx i, Idx j) {
	Idx it = (j * m_ni + i) * 4 + 2;
	Idx it0 = (j == m_nj - 1) ? static_cast<Idx>(Boundaries::SOLID_WALL) :
		4 * ((j + 1) * m_ni + i) + 0;
	Idx it1 = (j * m_ni + i) * 4 + 3;
	Idx it2 = (j * m_ni + i) * 4 + 1;
	m_triang_triangs.row(it) << it0, it1, it2;
	Idx ip0 = (m_ni + 1) * (j + 1) + i + 1;
	Idx ip1 = (m_ni + 1) * (j + 1) + i;
	Idx ip2 = m_num_verts + m_ni * j + i;
	m_triang_points.row(it) << ip0, ip1, ip2;
	m_edge_points.row(m_curr_e) << ip1, ip0;
	m_edge_triangs.row(m_curr_e++) << it, it0;
	m_edge_points.row(m_curr_e) << ip1, ip2;
	m_edge_triangs.row(m_curr_e++) << it, it1;
	m_triang_edges.row(it) << m_curr_e - 2, m_curr_e - 1, m_curr_e - 3;
}

void StructTriangMesh::TriangL(Idx i, Idx j) {
  Idx it = (j * m_ni + i) * 4 + 3;
  Idx it0 = (i == 0) ? static_cast<Idx>(Boundaries::SOLID_WALL) :
    4 * (j * m_ni + i - 1) + 1;
  Idx it1 = (j * m_ni + i) * 4 + 0;
  Idx it2 = (j * m_ni + i) * 4 + 2;
  m_triang_triangs.row(it) << it0, it1, it2;
  Idx ip0 = (m_ni + 1) * (j + 1) + i;
  Idx ip1 = (m_ni + 1) * j + i;
  Idx ip2 = m_num_verts + m_ni * j + i;
  m_triang_points.row(it) << ip0, ip1, ip2;
  if (i == 0) {
    // we add left edge of this square-box
    m_edge_points.row(m_curr_e) << ip1, ip0;
    m_edge_triangs.row(m_curr_e++) << it, it0;
    // so now our left edge is exactly
    // m_curr_e-1
    m_triang_edges.row(it) << m_curr_e - 1, m_curr_e - 6, m_curr_e - 2;
    return;
  }
  if (i == 1) {
    // but in this case we already have edge
    // with this coordinates
    // it's right edge of square {0, j}
    // so we should decrement m_curr_e in other way
    if (j == 0)
      m_triang_edges.row(it) << m_curr_e - 12, m_curr_e - 5, m_curr_e - 1;
    else
      m_triang_edges.row(it) << m_curr_e - 11, m_curr_e - 5, m_curr_e - 1;
    return;
  }
  // in general we do same thing as
  // in case (i == 1), but as the 
  // previous square had number {i-1, j}
  // and i-1 != 0, we now should decrement
  // m_curr_e on smaller value then in case
  // (i == 1)
  if (j == 0)
    m_triang_edges.row(it) << m_curr_e - 11, m_curr_e - 5, m_curr_e - 1;
  else
    m_triang_edges.row(it) << m_curr_e - 10, m_curr_e - 5, m_curr_e - 1;
}

