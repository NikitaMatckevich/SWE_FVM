#pragma once
#include "ValueField.h"

struct MUSCLObject {
  Bathymetry const& B_;
  VolumeField& Vol_;
  EdgeField& Edg_;
  
  using NodeArray = Eigen::ArrayXd;
  NodeArray max_wp_;

  void update();
private:
  bool dry_cell(index i) const;
  bool full_wet_cell(index i) const;
  bool part_wet_cell(index i) const;
  Eigen::Matrix32d gradient(index i) const;

  void full_wet_reconst(index i);
};