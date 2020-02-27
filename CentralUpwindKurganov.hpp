///////////////////////////////////////////////////////////
///                                                     ///
///   SECOND-ORDER FINITE-VOLUME SEMI-DISCRETE SOLVER   ///
///   FOR SHALLOW WATER FLOWS WITH WETTING AND DRYING   ///
///          OVER IRREGULAR BOTTOM TOPOGRAPHIES         ///
///                                                     ///
///           ICT SB RAS, Novosibirsk, Russia           ///
///                                                     ///
///////////////////////////////////////////////////////////

#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array> //TODO: replace all std containers with Eigen
#include "ConfigParser.hpp"
#include "DimensionManager.hpp"
#include "CubicPolyMath.hpp"
#include "TriangMesh.hpp"
#include "TriangSurfReconst.hpp"

using namespace Eigen;
using var2 = std::array<Array3d, 2>; //TODO: Eigen Tensor
namespace {
  bool is_wet(double h) {
    const double h_e = 1e-13; // minimal posible water depth of "wet" cell  
    return h > h_e;
  }
  //If cell isn't wet, velocities and depth should be equal to 0 
  Array3d dry_state(double B) { return { B, 0., 0. }; }
  //In MUSCL reconstruction we use 2 types of variables
  //Array3d Y is a vector with 3 "primitive" flow components:
  //  1) free suface elevation "w" or some f(w, R, phi, teta);
  //  2) 2D velocity vector    "u" or some f(u, R, phi, teta);
  //Array3d U is a similar vector in "conservative" form:
  //  1) free suface elevation "w" or some f(w, R, phi, teta);
  //  2) 2D impulse vector     "p" or some f(p, R, phi, teta);
  //Whenever it doesn't matter whether U or Y was passed as a
  //function argument the variable V is used (both U-type and
  //Y-type vectors may appear on the place of V)
  //B is bathymetry level
  //h is water depth: h = w - B
  inline double h(double V0, double B, double cos = 1.) {
    return V0 - B;
  }
  inline double w(double V0, double B = 0., double cos = 1.) {
    return V0;
  }
  inline Vector2d u(Array3d const & Y, double cos = 1.) {
    return Y.tail(2).matrix();
  }
  //Solver variables
  inline Array3d U(double w, double B, Vector2d u, double cos = 1.) {
    double h = w - B;
    Array3d res;
    if (is_wet(h)) {
      res[0] = w;
      res.tail(2) = h * u.array();
    }
    else res = dry_state(B);
    return res;
  }
  inline Array3d U(Array3d const & Y, double B) {
    Array3d res;
    double hi = h(Y[0], B);
    if (is_wet(hi)) {
      res = Y;
      res.tail(2) *= hi;
    }
    else res = dry_state(B);
    return res;
  }
  inline Array3d Y(double w, double B, Vector2d u, double cos = 1.) {
    double h = w - B;
    Array3d res;
    if (is_wet(h)) {
      res[0] = w;
      res.tail(2) = u.array();
    }
    else res = dry_state(B);
    return res;
  }
  inline Array3d Y(Array3d const & U, double B) {
    const double eps = 1e-4;
    Array3d res = U;
    double hi = h(U[0], B);
    if (hi < eps) {
      double pow_1 = hi*hi*hi*hi;
      const double pow_2 = 1e-8;
      double ih = sqrt(2)*hi / sqrt(pow_1 + pow_2);
      res.tail(2) *= ih;
    }
    else {
      res.tail(2) /= hi;
    }
    return res;
  }
  //Local wavespeeds
  double a_in (Vector2d const & n, double B,
    Array3d const & Y1, Array3d const & Y2) {
    double c1 = sqrt(h(Y1[0], B)); double c2 = sqrt(h(Y2[0], B));
    double lam_in  = u(Y1).dot(n) - c1;
    double lam_out = u(Y2).dot(n) - c2;
    return -std::min<double>({ lam_in, lam_out, 0. });
  }
  double a_out(Vector2d const & n, double B,
    Array3d const & Y1, Array3d const & Y2) {
    double c1 = sqrt(h(Y1[0], B)); double c2 = sqrt(h(Y2[0], B));
    double lam_in  = u(Y1).dot(n) + c1;
    double lam_out = u(Y2).dot(n) + c2;
    return std::max<double>({ lam_in, lam_out, 0. });
  }
  //Flux
  Array3d Flux(Vector2d const & n, double B, Array3d const & Y) {
    Array3d res;
    double hn = h(Y[0], B);
    if (is_wet(hn)) {
      double un = u(Y).dot(n);
      res[0] = hn * un;
      res.tail(2) = (hn*un*u(Y) + n*hn*hn*0.5).array();
    }
    else res.setZero();
    return res;
  }
  //TVD limiter from Delis et al., 2011
  double LIM(double a, double b) {
    const double LIM_e = 1e-14; // sufficiently small constant from LIM limiter function
    if (std::signbit(a) != std::signbit(b)) return 0.;
    return ((a*a + LIM_e)*b + (b*b + LIM_e)*a) / (a*a + b*b + 2.*LIM_e);
  }
  //Count gradient coefficients (on compact 3-stencil)
  Matrix32d grad_coefs(Array23d const & p) {
    Matrix32d res;
    double d = 1. / det(p.col(2) - p.col(0), p.col(2) - p.col(1));
    res.col(0) = d * (crcShiftL<double>(p.row(1)) - crcShiftR<double>(p.row(1))).matrix();
    res.col(1) = d * (crcShiftR<double>(p.row(0)) - crcShiftL<double>(p.row(0))).matrix();
    return res;
  }
  std::tuple<Array3d, Array3d> TVD(
    Array3d const & Y1t, Array3d const & Y1e, Array3d const & Y1c,
    Array3d const & Y2t, Array3d const & Y2e, Array3d const & Y2c,
    double portion1) {
    double portion2 = 1. - portion1;
    //Returns TVD-reconstructed values at edge center
    //See A.I. Delis et al., 2011
    //for details of this TVD reconstruction

    Array3d Central1 = portion1 * (Y2t - Y1t);
    Array3d AcrossEdge1 = Y1c - Y1t;
    Array3d AlongEdge1 = Y1e - Y1c;
    Array3d Upwind1 = 2.*AcrossEdge1 - Central1;
    Array3d Y1t_new = Y1t + AlongEdge1 +
      Upwind1.binaryExpr(Central1, [](double a, double b) {return LIM(a, b); });
    Array3d Central2 = portion2 * (Y1t - Y2t);
    Array3d AcrossEdge2 = Y2c - Y2t;
    Array3d AlongEdge2 = Y2e - Y2c;
    Array3d Upwind2 = 2.*AcrossEdge2 - Central2;
    Array3d Y2t_new = Y2t + AlongEdge2 +
      Upwind2.binaryExpr(Central2, [](double a, double b) {return LIM(a, b); });

    auto is_new_better = abs(Y2t_new - Y1t_new) < abs(Y2t - Y1t);
    return std::make_tuple(
      is_new_better.select(Y1t_new, Y1t),
      is_new_better.select(Y2t_new, Y2t)
    );
  }
};

using namespace std;
class Solver {
private:
     //continuous piecewise linear bathymetry reconstruction, size = M_.num_nodes()
  TriangSurfReconst const & B_;
     //triangular mesh
  TriangMesh const & M_;
     //maximal level of free surface elevation at mesh node
     //size = M_.num_nodes()
  ArrayXd max_wp_;
     //free surface derivatives at edge centers, size of each = M_.num_edges * 2 
  ArrayX2d dwdx_, dwdy_;
     //variables, size of each = M_.num_tiangles() * 3
  Array3Xd Y_, U_, U1_, U2_;
     //source terms
  double f_, tau_;
     //MUSCL reconstruction:
     //values at edge centers, size = M_.num_edges() * 2 * 3
  std::vector<var2> Ye_; //TODO: Eigen Tensor
     //values at edge intersections with centroid-crossing lines,
     //size = M_.num_edges() * 2 * 3
  std::vector<var2> Yc_; //TODO: Eigen Tensor
     //Kurganov fluxes et edge centers
     //size = M_.num_edges() * 3
  Array3Xd F_;
     //Wave local characteristics
  const double CFL_;
  const double sqrt_max_area_;
  const double min_wavespeed_to_capture_;
  double       min_length_to_wavespeed_;

  bool      dry_cell(node_tag i) const {
    return !is_wet(h(i));
  }
  bool full_wet_cell(node_tag i) const {
    double max_Bp = M_.triang_points(i).unaryExpr([this](node_tag i) {
      return this->B_.p(i);
    }).maxCoeff();
    bool has_dry_neighbour = M_.triang_triangs(i).unaryExpr([this](node_tag i) {
      return (i >= 0) ? dry_cell(i) : false ;
    }).any();
    return (max_Bp >= w(i)) && !has_dry_neighbour;
  }
  bool part_wet_cell(node_tag i) const {
    return !dry_cell(i) && !full_wet_cell(i);
  }

  //Returns 1st order gradient approximation at the triangle midpoint
  Matrix32d grad_Y(node_tag i) const {
    assert(full_wet_cell(i));
    triang_tag const & it = M_.triang_triangs(i);
    triang_tag const & ie = M_.triang_edges(i);
    //first step of reconstruction: 3-stencil plane's gradient 
    Array23d points; // points needed to build reconstruction plane 
    Matrix3d values; // values at these points
    for (int k = 0; k < 3; k++) {
      if (it[k] == SOLID_WALL) {
        points.col(k) = M_.t(i);
        values.col(k) = Y_.col(i).matrix();
      } else
      if (it[k] == FREE_FLOW) {
        points.col(k) = M_.t(i);
        values.col(k) = Y_.col(i).matrix();
      }
      else
      if (full_wet_cell(it[k])) {
        points.col(k) = M_.t(it[k]);
        values.col(k) = Y_.col(it[k]).matrix();
      } else {
        points.col(k) = M_.e(ie[k]);
        if (i == M_.edge_triangs(ie[k])[0]) {
          values.col(k) = 0.5*(Y_.col(i) + Ye_[ie[k]][1]).matrix();
        }
        else {
          values.col(k) = 0.5*(Y_.col(i) + Ye_[ie[k]][0]).matrix();
        }
      }
    }
    Matrix32d res = values * grad_coefs(points);
    //Positivity preserving reconstruction. Returns values at triangle vertices
    triang_tag const & ip = M_.triang_points(i);
    auto Yx3 = Y_.col(i).replicate<1, 3>().matrix();
    auto tx3 = M_.t(i).replicate<1, 3>();
    Matrix3d Yp = Yx3 + res * (M_.p(ip) - tx3).matrix();
    Array3d Bp = B_.p(ip);
    Array3d hp = Yp.row(0).transpose().array() - Bp;
    auto are_wet  = hp.unaryExpr([](double h) { return is_wet(h); });
    auto both_dry = [&, this, i](int k, int l) mutable {
      Yp.col(l) = dry_state(Bp[l]);
      int m = 3 - k - l;
      Yp.col(m) = 3. * this->Y_.col(i);
      Yp(0, m) -= Bp[k] + Bp[l];
      return Yp * grad_coefs(M_.p(ip));
    };
    if (!(are_wet.all())) {
      for (int k = 0; k < 3; k++) {
        if (!are_wet[k]) {
          Yp.col(k) = dry_state(Bp[k]);
          int k1 = (k + 1) % 3;
          int k2 = (k + 2) % 3;
          if (!are_wet[k1]) res = both_dry(k, k1); else
          if (!are_wet[k2]) res = both_dry(k, k2); else
          res *= 3.*h(i) / (3.*h(i) - hp[k]);
          break;
        }
      }
    }
    return res;
  }

  //Changes Ye_, max_wp_
  void part_wet_reconst_1order(node_tag i) {
    assert(part_wet_cell(i));
    triang_tag const & p = M_.triang_points(i);
    double B13 = B_.p(p).maxCoeff();
    double B23 = B_.p(p).minCoeff();
    double B12 = 3.*B_.t(i) - B23 - B13;
    double B_delimiter = B12 + 1./3.*(B13 - B12)*(B13 - B12) / (B13 - B23);
    double w_rec;
    if (w(i) <= B_delimiter) {
      w_rec = B23 + cbrt(3.*h(i)*(B13 - B23)*(B12 - B23));
    }
    else {
      double coef_2nd_degree = -3.*B13;
      double coef_1st_degree = 3.*(B12*B13 + B13*B23 - B12*B23);
      double coef_free = (B13 - B23) * (3.*h(i)*(B13 - B12) - B12*(B12 + B23))
                       - B23*B23*B13;
      w_rec = bisection(
        cubic::poly(coef_2nd_degree, coef_1st_degree, coef_free),
        /*find from*/ B12, /*to*/ B13);
    }
    triang_tag const & e = M_.triang_edges(i);
    for (int k = 0; k < 3; k++) {
      Array3d Yek = ::Y(w_rec, B_.e(e[k]), u(i));
      max_wp_[p[k]] = std::max(w_rec, max_wp_[p[k]]);
      if (i == M_.edge_triangs(e[k])[0]) {
        Ye_[e[k]][0] = Yek;
        dwdx_(e[k], 0) = 0.;
        dwdy_(e[k], 0) = 0.;
      }
      else {
        Ye_[e[k]][1] = Yek;
        dwdx_(e[k], 1) = 0.;
        dwdy_(e[k], 1) = 0.;
      }
    }
  }

  //Changes Ye_, Yc_, dwdx_, dwdy_, max_wp_
  void full_wet_reconst       (node_tag i) {
    assert(full_wet_cell(i));
    triang_tag const & p = M_.triang_points(i);
    triang_tag const & e = M_.triang_edges(i);
    triang_tag const & t = M_.triang_triangs(i);
    Array3d Ye, Yc;
    auto dY = grad_Y(i);
    for (int k = 0; k < 3; k++) {
      double wp;
      if (M_.is_edge_boundary(e[k])) {
        wp = ::w(Y_(0, i));
        if (t[k] == SOLID_WALL) {
          Ye_[e[k]][0] = ::Y(::w(Y_(0, i)), B_.e(e[k]), Vector2d::Zero());
          dwdx_(e[k], 0) = dY(0, 0);
          dwdy_(e[k], 0) = dY(0, 1);
        }
      }
      else {
        Ye = Y_.col(i) + (dY * (M_.e(e[k]) - M_.t(i)).matrix()).array();
        Yc = Y_.col(i) + (dY * (M_.c(e[k]) - M_.t(i)).matrix()).array();
        wp = ::w(Y_(0, i) + dY.row(0)*(M_.p(p[k]) - M_.t(i)).matrix());
        if (i == M_.edge_triangs(e[k])[0]) {
          Ye_[e[k]][0] = Ye;
          Yc_[e[k]][0] = Yc;
          dwdx_(e[k], 0) = dY(0, 0);
          dwdy_(e[k], 0) = dY(0, 1);
        }
        else {
          Ye_[e[k]][1] = Ye;
          Yc_[e[k]][1] = Yc;
          dwdx_(e[k], 1) = dY(0, 0);
          dwdy_(e[k], 1) = dY(0, 1);
        }
      }
      max_wp_[p[k]] = std::max(wp, max_wp_[p[k]]);
    }
  }

  //Changes Ye_, dwdx_, dwdy_, max_wp_
  void part_wet_reconst_2order(node_tag i) {
    assert(part_wet_cell(i));
    triang_tag p = M_.triang_points(i);
    //Sorting p by bathymetry elevation level
    if (B_.p(p[0]) > B_.p(p[1])) std::swap(p[0], p[1]);
    if (B_.p(p[1]) > B_.p(p[2])) std::swap(p[1], p[2]);
    if (B_.p(p[0]) > B_.p(p[1])) std::swap(p[0], p[1]);
    double w23 = max_wp_[p[0]];
    double h23 = w23 - B_.p(p[0]);
    double ratio_B = (B_.p(p[1]) - B_.p(p[0])) / (B_.p(p[2]) - B_.p(p[0]));
    double ratio_h = h(i) / h23;
    double h_delimiter1 = 1./3. * h23 * ratio_B;
    double h_delimiter2 = 1./3. * h23 * (2. - ratio_B);
    triang_tag e = M_.triang_edges(i);
    //Further, we will distinguish 3 kinds of edges:
    // 1) "lowest" edge e[0] = p[0], e[1] = p[1];
    // 2) "middle" edge e[0] = p[0], e[1] = p[2];
    // 3) "highest" one e[0] = p[1], e[1] = p[2];
    if ((M_.edge_points(e[0]) > M_.edge_points(e[1])).any()) std::swap(e[0], e[1]);
    if ((M_.edge_points(e[1]) > M_.edge_points(e[2])).any()) std::swap(e[1], e[2]);
    if ((M_.edge_points(e[0]) > M_.edge_points(e[1])).any()) std::swap(e[0], e[1]);
    Array33d Ye;
    Matrix<double, 1, 2> dwet;
    Matrix<double, 1, 3> grad_values;
    grad_values[0] = w23;
    Array23d grad_points;
    grad_points.col(0) = M_.p(p[0]);
    auto semiwet_edge = [this, i, &e, &Ye] (int k,
      double val1, double val2, double portion) {
      if (portion < 0.5) {
        Ye.col(k) = dry_state(B_.e(e[k]));
      }
      else {
        double c = 1. / (2. * portion);
        double w_rec = c * val1 + (1. - c) * val2;
        Ye.col(k) = ::Y(w_rec, B_.e(e[k]), u(i));
      }
    };
    //Case 1
    if (h(i) <= h_delimiter1) {
      Ye.col(2) = dry_state(B_.e(e[2]));
      double k2 = sqrt(3. * ratio_h / ratio_B);
      grad_values    [1] = k2 * B_.p(p[1]) + (1. - k2)*B_.p(p[0]);
      grad_points.col(1) = k2 * M_.p(p[1]) + (1. - k2)*M_.p(p[0]);
      semiwet_edge(0, grad_values[0], grad_values[1], k2);
      double k3 = sqrt(3. * ratio_h * ratio_B);
      grad_values    [2] = k3 * B_.p(p[2]) + (1. - k3)*B_.p(p[0]);
      grad_points.col(2) = k3 * M_.p(p[2]) + (1. - k3)*M_.p(p[0]);
      semiwet_edge(1, grad_values[0], grad_values[2], k3);
      dwet = grad_values * grad_coefs(grad_points);
    }
    //Case 3
    else if (h(i) >= h_delimiter2) {
      double delta_w = 1.5 * (h(i) - h_delimiter2);
      grad_values[1] = delta_w + B_.p(p[1]) + (1. - ratio_B)*h23;
      grad_points.col(1) = M_.p(p[1]);
      grad_values[2] = delta_w + B_.p(p[2]);
      grad_points.col(2) = M_.p(p[2]);
      dwet = grad_values * grad_coefs(grad_points);
      double w_rec;
      w_rec = 0.5 * (grad_values[0] + grad_values[1]);
      Ye.col(0) = ::Y(w_rec, B_.e(e[0]), u(i));
      w_rec = 0.5 * (grad_values[0] + grad_values[2]);
      Ye.col(1) = ::Y(w_rec, B_.e(e[1]), u(i));
      w_rec = 0.5 * (grad_values[1] + grad_values[2]);
      Ye.col(2) = ::Y(w_rec, B_.e(e[2]), u(i));
    }
    //Case 2
    else {
      double a = 3. * ratio_h;
      double b = 1. - ratio_B;
      double ib = 1. / b;
      double k1 = 1. -
        bisection(cubic::poly(ib*(a - 3.), ib*ib*(1. + b - a)));
      double k3 = 1. - b*(1. - k1);
      double Bp1 = k1*B_.p(p[2]) + (1. - k1)*B_.p(p[1]);
      double Bp3 = k3*B_.p(p[2]) + (1. - k3)*B_.p(p[0]);
      grad_values    [1] = B_.p(p[1]) + (k1 / k3) * b * h23;
      grad_points.col(1) = M_.p(p[1]);
      grad_values    [2] = Bp1;
      grad_points.col(2) = k1*M_.p(p[2]) + (1. - k1)*M_.p(p[1]);
      dwet = grad_values * grad_coefs(grad_points);
      double w_rec;
      w_rec = 0.5 * (grad_values[0] + grad_values[1]);
      Ye.col(0) = ::Y(w_rec, B_.e(e[0]), u(i));
      semiwet_edge(1, grad_values[0], Bp3, k3);
      semiwet_edge(2, grad_values[1], Bp1, k1);
    }
    for (int k = 0; k < 3; k++) {
      if (i == M_.edge_triangs(e[k])[0]) {
        Ye_[e[k]][0] = Ye.col(k);
        dwdx_(e[k], 0) = dwet[0];
        dwdy_(e[k], 0) = dwet[1];
      }
      else {
        Ye_[e[k]][1] = Ye.col(k);
        dwdx_(e[k], 1) = dwet[0];
        dwdy_(e[k], 1) = dwet[1];
      }
    }
  }

  //Changes max_wp_, Ye_, Yc_, dwdx_, dwdy_
  void    MUSCL_reconst() {
    max_wp_ = B_.buffer();
    size_t n = M_.num_triangles();
    for (node_tag i = 0; i < n; i++) {
      if (part_wet_cell(i)) part_wet_reconst_1order(i);
    }
    for (node_tag i = 0; i < n; i++) {
      if (full_wet_cell(i)) full_wet_reconst(i);
    }
    //for (node_tag i = 0; i < n; i++) {
    //  if (part_wet_cell(i)) part_wet_reconst_2order(i);
    //}
  }

  //Changes F_, min_length_to_wavespeed_
  void     Flux_reconst() {
    min_length_to_wavespeed_ = 1.;
    for (node_tag e = 0; e < M_.num_edges(); e++) {
      node_tag t0 = M_.edge_triangs(e)[0];
      node_tag t1 = M_.edge_triangs(e)[1];
      auto n = M_.normal(e);
      if (M_.is_edge_boundary(e)) {
        switch (t1) {
        case SOLID_WALL : F_.col(e) = Flux(n, B_.e(e), Ye_[e][0]);
        }
      }
      else if (dry_cell(t0) && dry_cell(t1)) {
        F_.col(e).setZero();
      }
      else {
        Array3d Y0, Y1;
        if (full_wet_cell(t0) && full_wet_cell(t1)) {
          double portion = length(M_.t(t0), M_.c(e)) / length(M_.t(t0), M_.t(t1));
          std::tie(Y0, Y1) = TVD(
            Y_.col(t0), Ye_[e][0], Yc_[e][0],
            Y_.col(t1), Ye_[e][1], Yc_[e][1],
            portion);
        }
        else {
          Y0 = Ye_[e][0]; Y1 = Ye_[e][1];
        }
        double ak_in  = a_in (n, B_.e(e), Y0, Y1);
        double ak_out = a_out(n, B_.e(e), Y0, Y1);
        if (ak_in + ak_out > min_wavespeed_to_capture_) {
          ///////////////////////////////////////////////////
          double r0 = 2. * M_.area(t0) / M_.l(e);
          double r1 = 2. * M_.area(t1) / M_.l(e);
          double length_to_wavespeed =
            std::min(r0, r1) / std::max(ak_in, ak_out);
          if (length_to_wavespeed < min_length_to_wavespeed_)
            min_length_to_wavespeed_ = length_to_wavespeed;
          ///////////////////////////////////////////////////
          F_.col(e) = (ak_in  * Flux(n, B_.e(e), Y1) +
                       ak_out * Flux(n, B_.e(e), Y0) +
                       ak_in*ak_out*(::U(Y0, B_.e(e)) - ::U(Y1, B_.e(e))))
                    / (ak_in + ak_out);
        }
        else
          F_.col(e).setZero();
      }
    }
  }

  //Changes Y_, U_
  void make_Y() {
    for (node_tag i = 0; i < M_.num_triangles(); i++) {
      if (dry_cell(i)) {
        Y_.col(i) = dry_state(B_.t(i));
        U_.col(i) = dry_state(B_.t(i));
      }
      else {
        Y_.col(i) = ::Y(U_.col(i), B_.t(i));
        U_.col(i).tail(2) = Y_.col(i).tail(2)*h(i);
      }
    }
  }

  double dt_k(double dt, node_tag i, node_tag k) const {
    auto dt_drain = [&](node_tag i)->double {
      triang_tag const & e = M_.triang_edges(i);
      double sum = 0.;
      for (int k = 0; k < 3; k++) {
        double f = (i == M_.edge_triangs(e[k])[0]) ?
           F_(0, e[k]) :
          -F_(0, e[k]) ;
        sum += std::max(0., f);
      }
      if (sum > 0.)  return (M_.area(i) * h(i) / sum);
      else           return INFINITY;
    };
    node_tag e = M_.triang_edges(i)[k];
    node_tag ik;
    double f;
    if (i == M_.edge_triangs(e)[0]) {
      ik =  M_.edge_triangs(e)[1];
      f  =  F_(0, e);
    }
    else {
      ik =  M_.edge_triangs(e)[0];
      f  = -F_(0, e);
    }
    if (f >= 0.) {
                   return std::min(dt, dt_drain(i));
    }
    else {
      if (ik >= 0) return std::min(dt, dt_drain(ik));
      else         return dt;
    }
  }
  Array3d RHS(node_tag i, double dt) const {
    Array3d res = { 0., f_*U_(2,i), -f_*U_(1,i) };
    for (int k = 0; k < 3; k++) {
      node_tag ek = M_.triang_edges (i)[k];
      node_tag pk = M_.triang_points(i)[k];
      double iS = 1. / M_.area(i);
      double l = M_.l(ek);
      double h, dwdx, dwdy;
      auto n = M_.normal(ek);
      Array3d F = F_.col(ek);
      if (i == M_.edge_triangs(ek)[0]) {
        h = ::h(Ye_[ek][0][0], B_.e(ek));
        dwdx = dwdx_(ek, 0);
        dwdy = dwdy_(ek, 0);
      }
      else {
        n *= -1.;
        F *= -1.;
        h = ::h(Ye_[ek][1][0], B_.e(ek));
        dwdx = dwdx_(ek, 1);
        dwdy = dwdy_(ek, 1);
      }
      double dtk = dt_k(dt, i, k);
      res    -= dtk * iS*l*F;
      res[1] += dtk * (0.5*iS*l*n[0]*h*h - 1./3.*dwdx*h);
      res[2] += dtk * (0.5*iS*l*n[1]*h*h - 1./3.*dwdy*h);
    }
    return res;
  }

  //Changes U_, U1_, U2_,
  //  max_wp_, dwdx_, dwdy_,
  //  Y_, Ye_, Yc_,
  //  F_,
  //  min_length_to_wavespeed_;
  void SSP_RK(double dt) {
    size_t nt = M_.num_triangles(); 
    double itau = 1. / (1. + tau_*dt); // implicit treatment of the tau term
    // R-K method
    double a21 = 0.75;   double a22 = 0.25;   double b2 = 0.25;
    double a31 = 1./3.;  double a32 = 2./3.;  double b3 = 2./3.;

    // 1)
    for (node_tag i = 0; i < nt; i++) {
      U1_.col(i) = U_.col(i) + RHS(i, dt);
    }
    U1_.bottomRows(2) *= itau;
    U_.swap(U1_);
    make_Y();
    MUSCL_reconst();
    Flux_reconst();

    // 2)
    for (node_tag i = 0; i < nt; i++) {
      U2_.col(i) = a21 * U1_.col(i) + a22 * U_.col(i) + RHS(i, b2*dt);
    }
    U2_.bottomRows(2) *= itau;
    U_.swap(U2_);
    make_Y();
    MUSCL_reconst();
    Flux_reconst();

    // 3)
    for (node_tag i = 0; i < nt; i++) {
      // we don't need 3rd stage vector U3_, as we can reuse U2_
      U2_.col(i) = a31 * U1_.col(i) + a32 * U_.col(i) + RHS(i, b3*dt);
    }
    U2_.bottomRows(2) *= itau;
    U_.swap(U2_);
    make_Y();
    MUSCL_reconst();
    Flux_reconst();
  }

public:
  //Constructor
  Solver(Parser       const & Par,
    DimensionManager  const & Dim,
    TriangSurfReconst const & Bat,
    const double CFL = 0.9,
    const double min_wavespeed_to_capture = 1e-8,
    double initial_min_length_to_wavespeed = 1.)
    : B_(Bat)
    , M_(B_.mesh())
    , max_wp_(B_.buffer())
    , f_(Dim.unscale<scales::source>(Par.get<double>("Common", "f")))
    , tau_(Dim.unscale<scales::source>(Par.get<double>("Common", "tau")))
    , CFL_(CFL)
    , sqrt_max_area_(sqrt(max_triang_area(M_)))
    , min_wavespeed_to_capture_(min_wavespeed_to_capture)
    , min_length_to_wavespeed_(initial_min_length_to_wavespeed) {
    try {
      if ((M_.num_nodes() < 3) || (M_.num_triangles() < 1)) {
        std::ostringstream message;
        message << "C-U K: Invalid Mesh in input";
        throw std::invalid_argument(message.str());
      }
      if  (M_.num_nodes() != B_.size()) {
        std::ostringstream message;
        message << "C-U K: Invalid Bp vector in input";
        throw std::invalid_argument(message.str());
      }
      size_t ne = M_.num_edges();
      dwdx_.resize(ne, NoChange);
      dwdy_.resize(ne, NoChange);
      Ye_.resize(ne);
      Yc_.resize(ne);
      F_.resize(NoChange, ne);
      size_t nt = M_.num_triangles();
      Y_.resize(NoChange, nt);
      U_ .resize(NoChange, nt);
      U1_.resize(NoChange, nt);
      U2_.resize(NoChange, nt);
    }
    catch (std::invalid_argument const & e) {
      std::cerr << e.what() << std::endl;
    }
  }

  //Input U, changes U_, U1_, U2_
  int initialize_mareogram(Array3Xd const & input) {
    if (U_.size() != input.size()) {
      std::cerr << "Incorrect \"U\" input size in Solver" << std::endl;
      return -1;
    }
    U_  = input;
    U1_ = input;
    U2_ = input;
    make_Y();
    MUSCL_reconst();
    Flux_reconst();
    return 0;
  }

  //Dump
  Array3d  const & U(node_tag i) const { return U_.col(i); }
  Array3Xd const & U(/* all  */) const { return U_; }
  double           h(node_tag i) const { return U_(0, i) - B_.t(i); }
  double           w(node_tag i) const { return U_(0, i); }
  ArrayXd          w(/* all  */) const { return U_.row(0); }
  Vector2d         p(node_tag i) const { return U_.col(i).tail(2).matrix(); }
  Array2Xd         p(/* all  */) const { return U_.bottomRows(2); }

  Array3d  const & Y(node_tag i) const { return Y_.col(i); }
  Array3Xd const & Y(/* all  */) const { return Y_; }
  Vector2d         u(node_tag i) const { return Y_.col(i).tail(2).matrix(); }
  Vector2d         u(/* all  */) const { return Y_.bottomRows(2); }

  double CFL_dt() const { return 1./3.*0.5*CFL_*min_length_to_wavespeed_; }

  //see SSP_RK
  void run(double dt) { SSP_RK(1./3.*dt); }
};