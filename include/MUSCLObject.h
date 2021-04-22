#pragma once
#include "ValueField.h"

struct MUSCLObject {

  inline VolumeField& GetVolField() noexcept { return m_vol; }
  inline const VolumeField& GetVolField() const noexcept { return m_vol; }

	inline EdgeField& GetEdgField() noexcept { return m_edg; }
	inline const EdgeField& GetEdgField() const noexcept { return m_edg; }

	//inline VolumeField& GetSrcField() noexcept { return m_s; }
	//inline const VolumeField& GetSrcField() const noexcept { return m_s; }
  inline EdgeField& GetSrcField() noexcept { return m_src; }
	inline const EdgeField& GetSrcField() const noexcept { return m_src; }

	inline const Bathymetry& Bath() const noexcept { return m_b; };
  inline const TriangMesh& Mesh() const noexcept { return m_b.Mesh(); }

	bool IsDryCell    (Idx i) const;
  bool IsFullWetCell(Idx i) const;
  bool IsPartWetCell(Idx i) const;

	void ReconstructAll();
	void DumpFields(const std::string& filename) const;

 protected:

  MUSCLObject(Bathymetry&& b, const VolumeField& v0);
  MUSCLObject(Bathymetry&& b, VolumeField&& v0);

  Bathymetry m_b;

  VolumeField m_vol;
  EdgeField   m_edg;
	EdgeField   m_src;
  //VolumeField m_s;
  Storage<1>  m_max_w;
  
  Eigen::Matrix32d Gradient(Idx i) const;

	void ReconstructDryCell(Idx i);
  void ReconstructFullWetCell(Idx i);
	void ReconstructPartWetCell1(Idx i);
  void ReconstructPartWetCell2(Idx i);

};
