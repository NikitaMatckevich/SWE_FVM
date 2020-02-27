#pragma once
#include <iomanip>
#include "ConfigParser.hpp"
#include "DimensionManager.hpp"
#include "StructTriangMesh.hpp"
#include "TriangSurfReconst.hpp"
#include "CentralUpwindKurganov.hpp"

template <int rows>
using A_t = Eigen::Array<double, rows, Eigen::Dynamic>;
using U_t = Eigen::Array3Xd;
using V_t = Eigen::ArrayXd;

template <int rows>
A_t<rows> avg(A_t<rows> const & v, int Ni) {
  assert(v.cols() % 4 == 0);
  int avg_size = v.cols() / 4;
  assert(avg_size % Ni == 0);
  int Nj = avg_size / Ni;
  assert(Ni % 2 == 0);
  Ni /= 2;
  assert(Nj % 2 == 0);
  Nj /= 2;
  A_t<rows> avg_v(rows, avg_size);
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      int t = 4 * (j*Ni + i);
      int n1 = 16 * j * Ni + 8 * i;
      int n2 = n1 + 8 * Ni;
      avg_v.col(t + 0) = 0.25 * (v.col(n1 + 0) + v.col(n1 + 1) + v.col(n1 + 4) + v.col(n1 + 7));
      avg_v.col(t + 1) = 0.25 * (v.col(n1 + 5) + v.col(n1 + 6) + v.col(n2 + 4) + v.col(n2 + 5));
      avg_v.col(t + 2) = 0.25 * (v.col(n2 + 1) + v.col(n2 + 2) + v.col(n2 + 6) + v.col(n2 + 7));
      avg_v.col(t + 3) = 0.25 * (v.col(n1 + 2) + v.col(n1 + 3) + v.col(n2 + 0) + v.col(n2 + 3));
    }
  }
  return avg_v;
}

double cmp_L1(V_t const & std, V_t const & num) {
  assert(std.size() == num.size());
  return abs(std - num).sum() / num.size();
}

int convergence_analyzer(Parser const & P, DimensionManager const & D) {
  const double t0 = 0.45;
  const int max_factor = 32;
  const int Ni0 = 5;
  const int Nj0 = 1;
  const double h0 = 0.2;
  auto sol = [=](int factor, U_t * U) {
    StructTriangMesh M(Ni0*factor, Nj0*factor, h0 / factor);
    V_t z(M.num_nodes());  z.fill(-1.);
    TriangSurfReconst B(M, z);
    if (factor == max_factor) {
      const double x0 = 0.5 * (M.min_x() + M.max_x());
      const double r = 0.1;
      const double A = 0.2;
      U->resize(Eigen::NoChange, M.num_triangles());
      U->bottomRows(2) = 0.;
      for (node_tag i = 0; i < M.num_triangles(); i++) {
        double x = M.t(i)[0];
        double a = (x - x0)*(x - x0);
        double w = A * exp(-a / (2.*r*r));
        U->operator()(0, i) = std::max(w, B.t(i));
      }
    }
    Solver L(P, D, B);
    L.initialize_mareogram(*U);
    const double dt = 0.4 * 1e-4;
    for (double t = 0; t < t0; t += dt) {
      std::cout << "CFL = " << L.CFL_dt() << std::endl;
      std::cout << "t   = " << t << std::endl;
      L.run(dt);
    }
    return L.w();
  };
  U_t U0;
  V_t num = sol(max_factor, &U0);
  auto std = avg<1>(num, Ni0*max_factor);
        U0 = avg<3>(U0 , Ni0*max_factor);
  std::ofstream fout("convergence.dat");
  double err_prev = NAN;
  double err_curr;
  for (int k = max_factor / 2; k >= 4; k /= 2) {
    num = sol(k, &U0);
    fout << k << " " << Ni0 * k << " " << Nj0 * k << " " << num.size() << " ";
    err_curr = cmp_L1(std, num);
    fout << std::scientific << std::setprecision(10) << err_curr << " ";
    double rate = log2(err_prev) - log2(err_curr);
    fout << std::scientific << std::setprecision(10) << rate << std::endl;
    err_prev = err_curr;
    std = avg<1>(std, Ni0 * k);
    U0 = avg<3>(U0, Ni0 * k);
  }
  fout.close();
  return 0;
}

int error_analyzer(Parser const & P, DimensionManager const & D) {
  const double t0 = 0.45;
  const int max_factor = 32;
  const int Ni0 = 5;
  const int Nj0 = 1;
  const double h0 = 0.2;
  auto sol = [=](int factor, U_t * U) {
    StructTriangMesh M(Ni0*factor, Nj0*factor, h0 / factor);
    V_t z(M.num_nodes());  z.fill(-1.);
    TriangSurfReconst B(M, z);
    if (factor == max_factor) {
      const double x0 = 0.5 * (M.min_x() + M.max_x());
      const double r = 0.1;
      const double A = 0.2;
      U->resize(Eigen::NoChange, M.num_triangles());
      U->bottomRows(2) = 0.;
      for (node_tag i = 0; i < M.num_triangles(); i++) {
        double x = M.t(i)[0];
        double a = (x - x0)*(x - x0);
        double w = A * exp(-a / (2.*r*r));
        U->operator()(0, i) = std::max(w, B.t(i));
      }
    }
    Solver L(P, D, B);
    L.initialize_mareogram(*U);
    const double dt = 0.4 * 1e-4;
    for (double t = 0; t < t0; t += dt) {
      std::cout << "CFL = " << L.CFL_dt() << std::endl;
      std::cout << "t   = " << t << std::endl;
      L.run(dt);
    }
    return L.w();
  };
  U_t U0;
  V_t num = sol(max_factor, &U0);
  auto std = avg<1>(num, Ni0*max_factor);
  U0 = avg<3>(U0, Ni0*max_factor);
  std::ofstream fout("convergence.dat");
  double err_prev = NAN;
  double err_curr;
  for (int k = max_factor / 2; k >= 4; k /= 2) {
    num = sol(k, &U0);
    fout << k << " " << Ni0 * k << " " << Nj0 * k << " " << num.size() << " ";
    err_curr = cmp_L1(std, num);
    fout << std::scientific << std::setprecision(10) << err_curr << " ";
    double rate = log2(err_prev) - log2(err_curr);
    fout << std::scientific << std::setprecision(10) << rate << std::endl;
    err_prev = err_curr;
    std = avg<1>(std, Ni0 * k);
    U0 = avg<3>(U0, Ni0 * k);
  }
  fout.close();
  return 0;
}