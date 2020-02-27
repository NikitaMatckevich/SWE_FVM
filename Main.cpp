#include <iostream>
#include <vector>
using namespace std;

#include "ConfigParser.hpp"
#include "DimensionManager.hpp"
#include "UnstructuredMesh.hpp"
#include "CentralUpwindKurganov.hpp"

#include <Eigen/Core>
std::string const msh_file = "D:/backup/Homework/gmsh/tutorial/rectangle.msh";
std::string const Cfg = "D:/backup/Homework/Git/tests/Test/Test/config.ini";
/*
int task_mesh() {
  Mesh M;
  if (M.read(msh_file) != 0) return -1;
  std::ofstream fout("my_mesh.txt");
  fout << "Nodes:" << std::endl;
  for (int k = 0; k < M.num_nodes(); k++)
    fout << k << ":   " << M.p(k).transpose() << std::endl;
  fout << std::endl << std::endl << "Triangles:" << std::endl;
  for (int k = 0; k < M.num_triangles(); k++) {
    fout << k << ":   " << M.t(k).transpose() << std::endl;
    fout << ' ' << ": l=" << M.area(k) << std::endl;
  }
  fout << std::endl << std::endl << "Edges:" << std::endl;
  for (int k = 0; k < M.num_edges(); k++) {
    fout << k << ": e=" << M.e(k).transpose() << std::endl;
    fout << ' ' << "  c=" << M.c(k).transpose() << std::endl;
    fout << ' ' << "  l=" << M.l(k) << std::endl;
    fout << ' ' << "  n=" << M.normal(k).transpose() << std::endl;
  }
  fout.close();
  return 0;
}*/
int task_gauss() {
  Parser Par(Cfg);
  DimensionManager Dim(Par);
  Mesh M;
  if (M.read(msh_file) != 0) return -1;
  Eigen::ArrayXd Bp(M.num_nodes());
  for (node_tag i = 0; i < M.num_nodes(); i++) {
    Bp[i] = -1.0;
  }
  CentralUpwingKurganov L(Par, Dim, M, Bp);
  std::vector<var> U(M.num_triangles());
  for (node_tag i = 0; i < M.num_triangles(); i++) {
    double x = M.t(i)[0];
    double y = M.t(i)[1];
    double x0 = 0.5;
    double y0 = 0.5;
    double r = 0.1;
    double A = 0.1;
    //double arg = (y - y0)*(y - y0);
    //double arg = (x - x0)*(x - x0);
    double arg = (x - x0)*(x - x0) + (y - y0)*(y - y0);
    double eta = A * exp(-arg / (2.0*r*r));
    if (eta < L.Bt(i)) eta = L.Bt(i);
    U[i] = { eta, 0.0, 0.0 };
  }
  if (L.initialize_mareogram(U) != 0) return -1;
  int rec = 0;
  int frame = 0;
  double t = 0;
  std::ofstream fout;
  for (int k = 0; k <= 1000; k++) {
    if (k >= rec) {
      std::cout << "rec = " << rec << ",  t = " << t << std::endl;
      double mass = 0.0;
      fout.open("eta_" + std::to_string(frame) + ".dat");
      for (node_tag i = 0; i < M.num_triangles(); i++) {
        double x = M.t(i)[0]; double y = M.t(i)[1];
        fout << x << ' ' << y << ' ' << L.w(i) << std::endl;
        mass += L.h(i) * M.area(i);
      }
      fout.close();
      std::cout << "   mass = " << mass << std::endl;
      rec += 100;
      frame++;
    }
    double dt = L.CFL_dt();
    t += dt;
    L.run(dt);
  }
  return 0;
}

int main() {
  //task_mesh();
  task_gauss();
  std::cout << "All tasks done. Exit program?" << std::endl;
  getchar();
	return 0;
}