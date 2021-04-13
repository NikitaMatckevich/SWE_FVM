#pragma once
#include <ConfigParser.h>
#include <DimensionManager.h>
#include <TriangMesh.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif // defined pi value

// abstract
class Test {
protected:
  double mid_x_, mid_y_; // location
  double f_, tau_; // sources
  Test(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : f_  (Dim.unscale<scales::source>(Par.get("Common", "f"  )))
    , tau_(Dim.unscale<scales::source>(Par.get("Common", "tau")))
    , mid_x_(0.5*(M.min_x() + M.max_x()))
    , mid_y_(0.5*(M.min_y() + M.max_y())) {}

public:
  virtual double B(double x, double y) const = 0;

  virtual inline double u(double x, double y, double t) const = 0;
  virtual inline double v(double x, double y, double t) const = 0;
  virtual inline double h(double x, double y, double t) const = 0;
  inline double w(double x, double y, double t) const { return h(x, y, t) + B(x, y); }

  bool is_wet(double x, double y, double t) const { return h(x, y, t) > 1e-10; }
};

// 1) abstract
class BowlTest : public Test {
protected:
  double delta_;
  BowlTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : Test(Par, Dim, M)
    , delta_(Par.get("Common", "delta")) {}

public:
  double B(double x, double y) const override {
    return delta_*((x - mid_x_)*(x - mid_x_) + (y - mid_y_)*(y - mid_y_) - 1.0);
  }
};

// 1.1) concrete
class BallTest : public BowlTest {
private:
  double x0_, y0_, u0_, v0_;
  double lam1_, lam2_, alpha_, beta_;
  double A1_, A2_, A3_, A4_;
public:
  BallTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : BowlTest(Par, Dim, M)
    , x0_(Dim.unscale<scales::length>(Par.get("Ball", "x0")) - mid_x_)
    , y0_(Dim.unscale<scales::length>(Par.get("Ball", "y0")) - mid_y_)
    , u0_(Dim.unscale<scales::velocity>(Par.get("Ball", "u0")))
    , v0_(Dim.unscale<scales::velocity>(Par.get("Ball", "v0")))
    , lam1_(2.  * delta_ - 0.25*tau_*tau_ + 0.25*f_*f_)
    , lam2_(0.5 * f_ * tau_) {
    alpha_ = lam1_ + sqrt(lam1_*lam1_ + lam2_*lam2_);
    beta_ = 0.5*sqrt(2.*alpha_);
    alpha_ = lam2_ / sqrt(2.*alpha_);
    A1_ = 0.5*(alpha_*alpha_ + beta_*beta_ + 0.5*alpha_*tau_ + 0.5*beta_*f_)*x0_
        + 0.25*(beta_*tau_ - alpha_*f_)*y0_
        + 0.5*(alpha_*u0_ + beta_*v0_);
    A1_ /= (alpha_ * alpha_ + beta_ * beta_);
    A2_ = 0.5*(alpha_*alpha_ + beta_*beta_ + 0.5*alpha_*tau_ + 0.5*beta_*f_)*y0_
        - 0.25*(beta_*tau_ - alpha_*f_)*x0_
        + 0.5*(alpha_*v0_ - beta_*u0_);
    A2_ /= (alpha_ * alpha_ + beta_ * beta_);
    A3_ = x0_ - A1_;
    A4_ = y0_ - A2_;
  }

  inline double phi_m(double t) const {
    return exp(-(0.5*tau_ - alpha_) * t) * (
      A1_ * cos((0.5*f_ - beta_) * t) + A2_ * sin((0.5*f_ - beta_) * t)
      );
  }
  inline double phi_p(double t) const  {
    return exp(-(0.5*tau_ + alpha_) * t) * (
      A3_ * cos((0.5*f_ + beta_) * t) + A4_ * sin((0.5*f_ + beta_) * t)
      );
  }
  inline double psi_m(double t) const {
    return exp(-(0.5*tau_ - alpha_) * t) * (
      A2_ * cos((0.5*f_ - beta_) * t) - A1_ * sin((0.5*f_ - beta_) * t)
      );
  }
  inline double psi_p(double t) const {
    return exp(-(0.5*tau_ + alpha_) * t) * (
      A4_ * cos((0.5*f_ + beta_) * t) - A3_ * sin((0.5*f_ + beta_) * t)
      );
  }
  inline double u(double x, double y, double t) const override {
    double res;
    is_wet(x, y, t) ?
      res = (0.5*f_ + beta_) * psi_p(t) - (0.5*tau_ + alpha_) * phi_p(t) +
            (0.5*f_ - beta_) * psi_m(t) - (0.5*tau_ - alpha_) * phi_m(t) :
      res = 0.0;
    return res;
  }
  inline double v(double x, double y, double t) const override {
    double res;
    is_wet(x, y, t) ?
      res =  - (0.5*f_ + beta_) * phi_p(t) - (0.5*tau_ + alpha_) * psi_p(t) -
               (0.5*f_ - beta_) * phi_m(t) - (0.5*tau_ - alpha_) * psi_m(t) :
      res = 0.0;
    return res;
  }
  inline double h(double x, double y, double t) const override {
  	double fx = x - mid_x_ - phi_p(t) - phi_m(t); 
		double fy = y - mid_y_ - psi_p(t) - psi_m(t);
	 	double res = delta_ * (1. - fx * fx - fy * fy);
    return std::max(0., res);
  }
};

// 1.2) abstract
class ThackerTest : public BowlTest {
protected:
  double H0, p0, q0;
  ThackerTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : BowlTest(Par, Dim, M)
    , H0(Dim.unscale<scales::height>(Par.get("Thacker", "H0")))
    , p0(Dim.unscale<scales::source>(Par.get("Thacker", "p0")))
    , q0(Dim.unscale<scales::source>(Par.get("Thacker", "q0"))) {}

public:
  virtual inline double p(double t) const = 0;
  virtual inline double q(double t) const = 0;
  double u(double x, double y, double t) const override {
    return p(t)*(x - mid_x_) + q(t)*(y - mid_y_);
  }
  virtual double v(double x, double y, double t) const = 0;
  virtual double Hc(double t) const = 0;
  virtual double Hxx(double t) const = 0;
  virtual double Hxy(double t) const = 0;
  virtual double Hyy(double t) const = 0;
  double h(double x, double y, double t) const override {
    double res = Hc(t) +
      0.5*Hxx(t)*(x - mid_x_)*(x - mid_x_) +
          Hxy(t)*(x - mid_x_)*(y - mid_y_) + 
      0.5*Hyy(t)*(y - mid_y_)*(y - mid_y_);
    return std::max(0., res);
  }
};

// 1.2.1) concrete
class ConservativeThackerTest : public ThackerTest {
private:
  double w;
  double A;
public:
  ConservativeThackerTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : ThackerTest(Par, Dim, M)
    , A(p0 * p0 + q0 * q0)
    , w(sqrt(0.25*f_*f_ + 2.0*A + 4.0*delta_) + 0.5*f_) {}

  inline double p(double t) const override {
    return p0 * cos(w*t) + q0 * sin(w*t);
  }
  inline double q(double t) const override {
    return q0 * cos(w*t) - p0 * sin(w*t);
  }
  double v(double x, double y, double t) const override {
    return q(t)*(x - mid_x_) + p(t)*(mid_y_ - y);
  }
  double Hc(double t) const override {
    return H0;
  }
  double Hxx(double t) const override {
    return -A - q(t)*(w - f_) - 2.0*delta_;
  }
  double Hxy(double t) const override {
    return p(t)*(w - f_);
  }
  double Hyy(double t) const override {
    return -A + q(t)*(w - f_) - 2.0*delta_;
  }
};

// 1.2.2) concrete
class SolenoidalThackerTest : public ThackerTest {
private:
  double w;
public:
  SolenoidalThackerTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : ThackerTest(Par, Dim, M)
    , w(sqrt(q0 - 2.0*delta_)) {
    if (f_ == 0) assert(p0 == q0 * (1.0 - delta_));
    if (tau_ == 0) assert((p0 == 0) || (p0 == sqrt(q0 - 2.0*delta_)));
  }
    
  inline double p(double t) const override {
    return p0 * exp(-tau_ * t);
  }
  inline double q(double t) const override {
    return q0 * exp(-tau_ * t);
  }
  double v(double x, double y, double t) const override {
    return q(t)*(mid_x_ - x) + p(t)*(mid_y_ - y);
  }
  double Hc(double t) const override {
    return H0;
  }
  double Hxx(double t) const override {
    if (tau_ != 0) return 0;
    else return q0 * q0 - p0 * p0 - q0 * f_ - 2.0*delta_;
  }
  double Hxy(double t) const override {
    return p0 * f_;
  }
  double Hyy(double t) const override {
    if (tau_ != 0) return 0;
    else return q0 * q0 - p0 * p0 - q0 * f_ - 2.0*delta_;
  }
};

// 1.2.3) concrete
class ClassicThackerTest : public ThackerTest {
private:
  double w;
  double A, B;
public:
  ClassicThackerTest(Parser const & Par, DimensionManager const & Dim, TriangMesh const & M)
    : ThackerTest(Par, Dim, M)
    , w(sqrt(f_*f_ + 8.0*delta_)) {
    double res;
    res = (q0 - 0.5*f_)*(q0 - 0.5*f_) + 2.0*H0*H0 + p0 * p0 - 0.25*w*w;
    A = sqrt(res*res + w * w*p0*p0) / (res + 0.5*w*w);
    B = atan(w*p0 / res);
  }

  inline double p(double t) const override {
    return 0.5*w * A*sin(w*t + B) / (1.0 - A * cos(w*t + B));
  }
  inline double q(double t) const override {
    return (q0 - 0.5*f_) * (1.0 - A * cos(B)) / (1.0 - A * cos(w*t + B)) + 0.5*f_;
  }
  double v(double x, double y, double t) const override {
    return q(t)*(mid_x_ - x) + p(t)*(y - mid_y_);
  }
  double Hc(double t) const override {
    return H0 * (1.0 - A * cos(B)) / (1.0 - A * cos(w*t + B));
  }
  double Hxx(double t) const override {
    return 0.25*(w*w*(A*A - 1.0) + 2.0*(q0 - f_)*(q0 - f_)*(1.0 - A * cos(B))*(1.0 - A * cos(B)))
      / (1.0 - A * cos(w*t + B));
  }
  double Hxy(double t) const override {
    return 0;
  }
  double Hyy(double t) const override {
    return 0.25*(w*w*(A*A - 1.0) + 2.0*(q0 - f_)*(q0 - f_)*(1.0 - A * cos(B))*(1.0 - A * cos(B)))
      / (1.0 - A * cos(w*t + B));
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
