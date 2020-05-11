#pragma once
#include "ValueField.h"

namespace utils {

  struct MUSCLObject {
    inline const VolumeField& currentVolField() const { return Vol_; }
    inline const EdgeField& currentEdgField() const { return Edg_; }
  protected:
    MUSCLObject(VolumeField&& v0);
    using NodeField = Eigen::ArrayXd;

    VolumeField Vol_;
    EdgeField   Edg_;
    NodeField   Max_W_;

    inline const Bathymetry& bathymetry() const { return Vol_.bathymetry(); };

    bool is_dry_cell(index i) const;
    bool is_full_wet_cell(index i) const;
    bool is_part_wet_cell(index i) const;
    Eigen::Matrix32d gradient(index i) const;

    void reconstruct_full_wet_cell(index i);
    void reconstruct_all();
  };
} // namespace utils 