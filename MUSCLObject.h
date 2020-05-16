#pragma once
#include "ValueField.h"

struct MUSCLObject {
  using Array = BaseValueField::Array;
  inline const VolumeField& currentVolField() const { return Vol_; }
  inline const EdgeField&   currentEdgField() const { return Edg_; }
protected:
  MUSCLObject(Bathymetry&& b, VolumeField&& v0);
  using NodeField = Eigen::ArrayXd;
  Bathymetry b_;
  VolumeField Vol_;
  EdgeField   Edg_;
  NodeField   Max_W_;

  inline const Bathymetry& bath() const { return b_; };
  inline const TriangMesh& mesh() const { return b_.mesh(); }

  bool is_dry_cell(index i) const;
  bool is_full_wet_cell(index i) const;
  bool is_part_wet_cell(index i) const;
  Eigen::Matrix32d gradient(index i) const;

  void reconstruct_full_wet_cell(index i);
  void reconstruct_all();
};