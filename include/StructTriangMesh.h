#pragma once
#include "TriangMesh.h" 

struct StructTriangMesh : public TriangMesh {
  StructTriangMesh(size_t ni, size_t nj, double h);
  inline size_t Ni() const { return m_ni; }
  inline size_t Nj() const { return m_nj; } 
 private:
  Idx m_curr_e = 0;
  size_t m_ni, m_nj, m_num_verts, m_num_sqrs;
  void TriangB(Idx i, Idx j);
	void TriangR(Idx i, Idx j);
	void TriangT(Idx i, Idx j);
  void TriangL(Idx i, Idx j);
};
