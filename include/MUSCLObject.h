#pragma once
#include "ValueField.h"

struct MUSCLObject {

  inline VolumeField& GetVolField() noexcept { return m_vol; }
  inline const VolumeField& GetVolField() const noexcept { return m_vol; }

	inline const Bathymetry& Bath() const noexcept { return m_b; };
  inline const TriangMesh& Mesh() const noexcept { return m_b.Mesh(); }

	bool IsDryCell    (Idx i) const;
  bool IsFullWetCell(Idx i) const;
  bool IsPartWetCell(Idx i) const;
	
 protected:

  MUSCLObject(Bathymetry&& b, const VolumeField& v0);
  MUSCLObject(Bathymetry&& b, VolumeField&& v0);

  struct MUSCL {
    MUSCL(const Bathymetry& b, const Array<3>& at_origin, const Eigen::Matrix32d& grad, Idx i)
      : m_b(b)
      , m_at_origin(at_origin)
      , m_grad(grad)
      , m_i(i)
      {}

    inline Idx OriginId() const noexcept { return m_i; }
    inline Array<3> AtOrigin() const noexcept { return m_at_origin; }
    inline Eigen::Matrix32d Gradient(const Point& p) const noexcept {
      double h = AtPoint(p)[0] - m_b.AtPoint(m_i, p);
      if (h >= 0) {
        return m_grad;
      } else {
        Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
        grad.row(0) = m_b.Gradient(m_i).matrix().transpose();
        return grad;
      }
    }
    inline Array<3> AtPoint(const Point& p) const noexcept {
      Array<3> at_point = m_at_origin + (m_grad * (p - m_b.Mesh().T(m_i)).matrix()).array();
      double b = m_b.AtPoint(m_i, p);
      if ((at_point[0] - b) >= 0)
        return at_point;
      else
        return DryState(b);
    }

   private:
    const Bathymetry& m_b;
    const Array<3> m_at_origin;
    const Eigen::Matrix32d m_grad;
    Idx m_i;
  };

  Bathymetry m_b;

  VolumeField m_vol;
  Storage<1>  m_max_wp;
  
	MUSCL ReconstructDryCell     (Idx i) const;
  MUSCL ReconstructFullWetCell (Idx i) const;
	MUSCL ReconstructPartWetCell1(Idx i) const;
  MUSCL ReconstructPartWetCell2(Idx i) const;
};
