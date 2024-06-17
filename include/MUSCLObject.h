#pragma once
#include "ValueField.h"

struct MUSCLObject {

    inline VolumeField& GetVolField() noexcept { return m_vol; }
    inline const VolumeField& GetVolField() const noexcept { return m_vol; }

    inline const Domain& GetDomain() const noexcept { return m_b; };

    bool IsDryCell    (Idx i) const;
    bool IsFullWetCell(Idx i) const;
    bool IsPartWetCell(Idx i) const;
	
    struct MUSCL {
        MUSCL(const Domain& b, const Array<3>& at_origin, const Eigen::Matrix32d& grad, Idx i)
            : m_b(b)
            , m_at_origin(at_origin)
            , m_grad(grad)
            , m_i(i)
        {}

        inline Idx OriginId() const noexcept { return m_i; }
        inline Array<3> AtOrigin() const noexcept { return m_at_origin; }
    
        inline Array<3> AtPoint(const Point& p) const noexcept {
            Array<3> at_point = m_at_origin + (m_grad * (p - m_b.T(m_i)).matrix().topRows(2)).array();
            if ((at_point[0] - p[2]) >= 0) {
                return at_point;
            } else {
                return DryState(p[2]);
            }
        }
    
        inline Eigen::Matrix32d dryGradient() const noexcept {
            Eigen::Matrix32d grad = Eigen::Matrix32d::Zero();
            grad.row(0) = m_b.TriangSlope(m_i).transpose();
            return grad;
        }
    
        inline Eigen::Matrix32d Gradient(const Point& p) const noexcept {
          double h = AtPoint(p)[0] - p[2];
          if (h >= 0) {
            return m_grad;
          } else {
            return dryGradient();
          }
        }

    private: 
        const Domain& m_b;
        const Array<3> m_at_origin;
        const Eigen::Matrix32d m_grad;
        Idx m_i;
    };

  MUSCL ReconstructDryCell     (Idx i) const;
  MUSCL ReconstructFullWetCell (Idx i) const;
  MUSCL ReconstructPartWetCell1(Idx i) const;
  MUSCL ReconstructPartWetCell2(Idx i) const;

protected:

  MUSCLObject(const Domain& b, const VolumeField& v0);

  const Domain& m_b;
  VolumeField m_vol;
  Storage<1>  m_max_wp;
};
