#pragma once
#include "ValueField.h"

struct MUSCLObject {

  inline VolumeField& getVolField() noexcept { return Vol_; }
  inline const VolumeField& getVolField() const noexcept { return Vol_; }
	inline EdgeField& getEdgField() noexcept { return Edg_; }
	inline const EdgeField&   getEdgField() const noexcept { return Edg_; }
	inline const Bathymetry& bath() const noexcept { return b_; };
  inline const TriangMesh& mesh() const noexcept { return b_.mesh(); }
	void reconstruct_all();
	void Dump(const std::string& filename) const;

 protected:
  MUSCLObject(Bathymetry&& b, const VolumeField& v0);
  MUSCLObject(Bathymetry&& b, VolumeField&& v0);

	using NodeField = Eigen::ArrayXd;
  Bathymetry b_;

  VolumeField Vol_;
  EdgeField   Edg_;
  NodeField   Max_W_;
  
  bool is_dry_cell     (idx i) const;
  bool is_full_wet_cell(idx i) const;
  bool is_part_wet_cell(idx i) const;
  Eigen::Matrix32d gradient(idx i) const;

  void reconstruct_full_wet_cell(idx i);
};
