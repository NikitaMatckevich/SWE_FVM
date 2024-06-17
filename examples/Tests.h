#pragma once
#include <TriangMesh.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif // defined pi value

// abstract
class Test {
protected:
  double m_mid_x, m_mid_y; // location
  double m_cor, m_tau; // sources
  Test(double mid_x, double mid_y)
    : m_cor(dim.Unscale<Scales::source>(par.Get("Common", "cor")))
    , m_tau(dim.Unscale<Scales::source>(par.Get("Common", "tau")))
    , m_mid_x(mid_x)
    , m_mid_y(mid_y) {}

public:
  virtual double b(double x, double y) const = 0;

  virtual inline double u(double x, double y, double t) const = 0;
  virtual inline double v(double x, double y, double t) const = 0;
  virtual inline double h(double x, double y, double t) const = 0;
  inline double w (double x, double y, double t) const { return h(x, y, t) + b(x, y); }
  inline double hu(double x, double y, double t) const { return h(x, y, t) * u(x, y, t); }
  inline double hv(double x, double y, double t) const { return h(x, y, t) * v(x, y, t); }

  bool IsWet(double x, double y, double t) const { return h(x, y, t) >= tol; }
};

class LakeAtRestTest : public Test {
public:
  LakeAtRestTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : Test(par, dim, mid_x, mid_y) {}

  inline double b(double x, double y) const override {
    return (1. < x) && (x < 3.) && (1. < y) && (y < 3.) ? -0.2 : -1.;
  }
  inline double u(double x, double y, double t) const { return 0.; }
  inline double v(double x, double y, double t) const { return 0.; }
  inline double h(double x, double y, double t) const { return std::max(0., -b(x,y)); }
};

// 1) abstract
class BowlTest : public Test {
protected:
  double m_delta;
  BowlTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : Test(par, dim, mid_x, mid_y)
    , m_delta(par.Get("Common", "delta")) {}

public:
  inline double b(double x, double y) const override {
    return m_delta*((x - m_mid_x)*(x - m_mid_x) + (y - m_mid_y)*(y - m_mid_y) - 1.0);
  }
};

// 1.1) concrete
class BallTest : public BowlTest {
private:
  double m_x0, m_y0, m_u0, m_v0;
  double m_lam1, m_lam2, m_alpha, m_beta;
  double m_a1, m_a2, m_a3, m_a4;
public:
  BallTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : BowlTest(par, dim, mid_x, mid_y)
    , m_x0(dim.Unscale<Scales::length>(par.Get("Ball", "x0")) - mid_x)
    , m_y0(dim.Unscale<Scales::length>(par.Get("Ball", "y0")) - mid_y)
    , m_u0(dim.Unscale<Scales::velocity>(par.Get("Ball", "u0")))
    , m_v0(dim.Unscale<Scales::velocity>(par.Get("Ball", "v0")))
    , m_lam1(2.  * m_delta - 0.25 * m_tau * m_tau + 0.25 * m_cor * m_cor)
    , m_lam2(0.5 * m_cor * m_tau)
  {
    m_alpha = m_lam1 + sqrt(m_lam1*m_lam1 + m_lam2*m_lam2);
    m_beta = 0.5*sqrt(2.*m_alpha);
    m_alpha = m_lam2 / sqrt(2.*m_alpha);
    m_a1 = 0.5*(m_alpha*m_alpha + m_beta*m_beta + 0.5*m_alpha*m_tau + 0.5*m_beta*m_cor)*m_x0
        + 0.25*(m_beta*m_tau - m_alpha*m_cor)*m_y0
        + 0.5*(m_alpha*m_u0 + m_beta*m_v0);
    m_a1 /= (m_alpha * m_alpha + m_beta * m_beta);
    m_a2 = 0.5*(m_alpha*m_alpha + m_beta*m_beta + 0.5*m_alpha*m_tau + 0.5*m_beta*m_cor)*m_y0
        - 0.25*(m_beta*m_tau - m_alpha*m_cor)*m_x0
        + 0.5*(m_alpha*m_v0 - m_beta*m_u0);
    m_a2 /= (m_alpha * m_alpha + m_beta * m_beta);
    m_a3 = m_x0 - m_a1;
    m_a4 = m_y0 - m_a2;
  }

  inline double PhiM(double t) const {
    return exp(-(0.5*m_tau - m_alpha) * t) * (
      m_a1 * cos((0.5*m_cor - m_beta) * t) + m_a2 * sin((0.5*m_cor - m_beta) * t)
      );
  }
  inline double PhiP(double t) const  {
    return exp(-(0.5*m_tau + m_alpha) * t) * (
      m_a3 * cos((0.5*m_cor + m_beta) * t) + m_a4 * sin((0.5*m_cor + m_beta) * t)
      );
  }
  inline double PsiM(double t) const {
    return exp(-(0.5*m_tau - m_alpha) * t) * (
      m_a2 * cos((0.5*m_cor - m_beta) * t) - m_a1 * sin((0.5*m_cor - m_beta) * t)
      );
  }
  inline double PsiP(double t) const {
    return exp(-(0.5*m_tau + m_alpha) * t) * (
      m_a4 * cos((0.5*m_cor + m_beta) * t) - m_a3 * sin((0.5*m_cor + m_beta) * t)
      );
  }
  inline double u(double x, double y, double t) const override {
    double res;
    IsWet(x, y, t) ?
      res = (0.5*m_cor + m_beta) * PsiP(t) - (0.5*m_tau + m_alpha) * PhiP(t) +
            (0.5*m_cor - m_beta) * PsiM(t) - (0.5*m_tau - m_alpha) * PhiM(t) :
      res = 0.0;
    return res;
  }
  inline double v(double x, double y, double t) const override {
    double res;
    IsWet(x, y, t) ?
      res =  - (0.5*m_cor + m_beta) * PhiP(t) - (0.5*m_tau + m_alpha) * PsiP(t) -
               (0.5*m_cor - m_beta) * PhiM(t) - (0.5*m_tau - m_alpha) * PsiM(t) :
      res = 0.0;
    return res;
  }
  inline double h(double x, double y, double t) const override {
  	double fx = x - m_mid_x - PhiP(t) - PhiM(t); 
		double fy = y - m_mid_y - PsiP(t) - PsiM(t);
	 	double res = m_delta * (1. - fx * fx - fy * fy);
    return std::max(0., res);
  }
};

// 1.2) abstract
class ThackerTest : public BowlTest {
protected:
  double H0, p0, q0;
  ThackerTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : BowlTest(par, dim, mid_x, mid_y)
    , H0(dim.Unscale<Scales::height>(par.Get("Thacker", "H0")))
    , p0(dim.Unscale<Scales::source>(par.Get("Thacker", "p0")))
    , q0(dim.Unscale<Scales::source>(par.Get("Thacker", "q0"))) {}

public:
  virtual inline double p(double t) const = 0;
  virtual inline double q(double t) const = 0;
  double u(double x, double y, double t) const override {
    return p(t)*(x - m_mid_x) + q(t)*(y - m_mid_y);
  }
  virtual double v(double x, double y, double t) const = 0;
  virtual double Hc(double t) const = 0;
  virtual double Hxx(double t) const = 0;
  virtual double Hxy(double t) const = 0;
  virtual double Hyy(double t) const = 0;
  double h(double x, double y, double t) const override {
    double res = Hc(t) +
      0.5*Hxx(t)*(x - m_mid_x)*(x - m_mid_x) +
          Hxy(t)*(x - m_mid_x)*(y - m_mid_y) + 
      0.5*Hyy(t)*(y - m_mid_y)*(y - m_mid_y);
    return std::max(0., res);
  }
};

// 1.2.1) concrete
class ConservativeThackerTest : public ThackerTest {
private:
  double w;
  double a;
public:
  ConservativeThackerTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : ThackerTest(par, dim, mid_x, mid_y)
    , a(p0 * p0 + q0 * q0)
    , w(sqrt(0.25*m_cor*m_cor + 2.0*a + 4.0*m_delta) + 0.5*m_cor) {}

  inline double p(double t) const override {
    return p0 * cos(w*t) + q0 * sin(w*t);
  }
  inline double q(double t) const override {
    return q0 * cos(w*t) - p0 * sin(w*t);
  }
  double v(double x, double y, double t) const override {
    return q(t)*(x - m_mid_x) + p(t)*(m_mid_y - y);
  }
  double Hc(double t) const override {
    return H0;
  }
  double Hxx(double t) const override {
    return -a - q(t)*(w - m_cor) - 2.0*m_delta;
  }
  double Hxy(double t) const override {
    return p(t)*(w - m_cor);
  }
  double Hyy(double t) const override {
    return -a + q(t)*(w - m_cor) - 2.0*m_delta;
  }
};

// 1.2.2) concrete
class SolenoidalThackerTest : public ThackerTest {
private:
  double w;
public:
  SolenoidalThackerTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : ThackerTest(par, dim, mid_x, mid_y)
    , w(sqrt(q0 - 2.0*m_delta))
  {
    if (m_cor == 0) assert(p0 == q0 * (1.0 - m_delta));
    if (m_tau == 0) assert((p0 == 0) || (p0 == sqrt(q0 - 2.0*m_delta)));
  }
    
  inline double p(double t) const override {
    return p0 * exp(-m_tau * t);
  }
  inline double q(double t) const override {
    return q0 * exp(-m_tau * t);
  }
  double v(double x, double y, double t) const override {
    return q(t)*(m_mid_x - x) + p(t)*(m_mid_y - y);
  }
  double Hc(double t) const override {
    return H0;
  }
  double Hxx(double t) const override {
    if (m_tau != 0) return 0;
    else return q0 * q0 - p0 * p0 - q0 * m_cor - 2.0*m_delta;
  }
  double Hxy(double t) const override {
    return p0 * m_cor;
  }
  double Hyy(double t) const override {
    if (m_tau != 0) return 0;
    else return q0 * q0 - p0 * p0 - q0 * m_cor - 2.0*m_delta;
  }
};

// 1.2.3) concrete
class ClassicThackerTest : public ThackerTest {
private:
  double w;
  double a, b;
public:
  ClassicThackerTest(const Parser& par, const DimensionManager& dim, double mid_x, double mid_y)
    : ThackerTest(par, dim, mid_x, mid_y)
    , w(sqrt(m_cor * m_cor + 8. * m_delta))
  {
    double qz = (q0 - 0.5 * m_cor) * (q0 - 0.5 * m_cor);
    double rz = qz + 2. * H0 * H0 + p0 * p0 - 0.25 * w * w;
    a = sqrt(rz * rz + w * w * p0 * p0) / (rz + 0.5 * w * w);
    b = atan(w * p0 / rz);

    double cp = (2. * q0 - m_cor) * (1. - a * cos(b)) / w;
    double cw = H0 * (1. - a * cos(b)); 
    std::cout << "cp=" << cp << ",  cw=" << cw << std::endl;
  }

  inline double p(double t) const override {
    return 0.5 * w * a * sin(w*t + b) / (1. - a * cos(w*t + b));
  }
  inline double q(double t) const override {
    return (q0 - 0.5 * m_cor) * (1. - a * cos(b)) / (1. - a * cos(w*t + b)) + 0.5 * m_cor;
  }
  double v(double x, double y, double t) const override {
    return q(t)*(m_mid_x - x) + p(t)*(y - m_mid_y);
  }
  double Hc(double t) const override {
    return H0 * (1. - a * cos(b)) / (1. - a * cos(w*t + b));
  }
  double Hxx(double t) const override {
    double qz = (q0 - 0.5 * m_cor) * (q0 - 0.5 * m_cor);
    double az0 = (1. - a * cos(b)) * (1. - a * cos(b));
    double azt = (1. - a * cos(w*t + b)) * (1. - a * cos(w*t + b));
    return (0.25 * w * w * (a * a - 1.) + qz * az0) / azt;
  }
  double Hxy(double t) const override {
    return 0;
  }
  double Hyy(double t) const override {
    return Hxx(t);
  }
};

/*
double u_triang(Eigen::Array23d const & triang, double t) const override {
  double S = triangle_area(triang.col(0), triang.col(1), triang.col(2));
  double x1 = triang(0, 0) - mid_x_;
  double y1 = triang(1, 0) - mid_y_;
  double x2 = triang(0, 1) - mid_x_;
  double y2 = triang(1, 1) - mid_y_;
  double x3 = triang(0, 2) - mid_x_;
  double y3 = triang(1, 2) - mid_y_;
  double Ix = x1 + x2 + x3;
  double Iy = y1 + y2 + y3;
  return (Ix*p(t) + Iy * q(t)) * S / 3.;
}
double v_triang(Eigen::Array23d const & triang, double t) const override {
  double S = triangle_area(triang.col(0), triang.col(1), triang.col(2));
  double x1 = triang(0, 0) - mid_x_;
  double y1 = triang(1, 0) - mid_y_;
  double x2 = triang(0, 1) - mid_x_;
  double y2 = triang(1, 1) - mid_y_;
  double x3 = triang(0, 2) - mid_x_;
  double y3 = triang(1, 2) - mid_y_;
  double Ix = x1 + x2 + x3;
  double Iy = y1 + y2 + y3;
  return (-Ix * q(t) + Iy * p(t)) * S / 3.;
}
double H_triang(Eigen::Array23d const & triang, double t) const override {
  double x1 = triang(0, 0) - mid_x_;
  double y1 = triang(1, 0) - mid_y_;
  double x2 = triang(0, 1) - mid_x_;
  double y2 = triang(1, 1) - mid_y_;
  double x3 = triang(0, 2) - mid_x_;
  double y3 = triang(1, 2) - mid_y_;
  double Ixx = (x1*x1 + x2 * x2 + x3 * x3 + x1 * x2 + x1 * x3 + x2 * x3);
  double Ixy = (x1*(2.*y1 + y2 + y3) +
    x2 * (2.*y2 + y1 + y3) +
    x3 * (2.*y3 + y1 + y2));
  double Iyy = (y1*y1 + y2 * y2 + y3 * y3 + y1 * y2 + y1 * y3 + y2 * y3);
  return (12.*Hc(t) + Ixx * Hxx(t) + Ixy * Hxy(t) + Iyy * Hyy(t)) / 12.;
}
*/
