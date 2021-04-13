#include <DimensionManager.h>
#include <SpaceIntegrators.h>
#include <TimeIntegrators.h>
#include <iostream>
#include <Tests.hpp>
#include <omp.h>

struct StructTriangMesh : public TriangMesh {
  StructTriangMesh(size_t Ni, size_t Nj, double h)
    : Ni_(Ni), Nj_(Nj)
    , num_verts_((Ni + 1)* (Nj + 1))
    , num_sqrs_(Ni* Nj) {
    size_t num_triangs = 4 * num_sqrs_;
    triang_points_.resize(num_triangs, Eigen::NoChange);
    triang_edges_.resize(num_triangs, Eigen::NoChange);
    triang_triangs_.resize(num_triangs, Eigen::NoChange);
    size_t num_pts = num_verts_ + num_sqrs_;
    p_.resize(Eigen::NoChange, num_pts);
    size_t num_edges = num_triangs + num_pts - 1;
    edge_points_.resize(num_edges, Eigen::NoChange);
    edge_triangs_.resize(num_edges, Eigen::NoChange);
    // fill p_ array with square's vertices
    for (idx j = 0; j < Nj_ + 1; j++) {
      for (idx i = 0; i < Ni_ + 1; i++) {
        p_(0, (Ni_ + 1) * j + i) = i * h;
        p_(1, (Ni_ + 1) * j + i) = j * h;
      }
    }
    // fill p_ array with square's centers
    for (idx j = 0; j < Nj_; j++) {
      for (idx i = 0; i < Ni_; i++) {
        p_(0, num_verts_ + Ni_ * j + i) = (i + 0.5) * h;
        p_(1, num_verts_ + Ni_ * j + i) = (j + 0.5) * h;
      }
    }
    // determine triangles and edges
    for (idx j = 0; j < Nj_; j++) {
      for (idx i = 0; i < Ni_; i++) {
        triang_b(i, j); //-bottom triangle of this square
        triang_r(i, j); //-right  triangle of this square
        triang_t(i, j); //-top    triangle of this square
        triang_l(i, j); //-left   triangle of this square
      }
    }
  }
  size_t Ni() const { return Ni_; }
  size_t Nj() const { return Nj_; }
 
 private:
  idx curr_e_ = 0;
  size_t Ni_, Nj_, num_verts_, num_sqrs_;
  void triang_b(idx i, idx j) {
    idx it = (j * Ni_ + i) * 4 + 0;
    idx it0 = (j == 0) ? static_cast<idx>(boundaries::SOLID_WALL) :
      4 * ((j - 1) * Ni_ + i) + 2;
    idx it1 = (j * Ni_ + i) * 4 + 1;
    idx it2 = (j * Ni_ + i) * 4 + 3;
    triang_triangs_.row(it) << it0, it1, it2;
    idx ip0 = (Ni_ + 1) * j + i;
    idx ip1 = (Ni_ + 1) * j + i + 1;
    idx ip2 = num_verts_ + Ni_ * j + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    if (j == 0) {
      edge_points_.row(curr_e_) << ip0, ip1;
      edge_triangs_.row(curr_e_++) << it, it0;
    }
    edge_points_.row(curr_e_) << ip1, ip2;
    edge_triangs_.row(curr_e_++) << it, it1;
    edge_points_.row(curr_e_) << ip0, ip2;
    edge_triangs_.row(curr_e_++) << it, it2;
    if (j == 0) {
      // we've already put all 3 edges in mesh
      triang_edges_.row(it) << curr_e_ - 3, curr_e_ - 2, curr_e_ - 1;
      return;
    }
    if (j == 1) {
      // take bottom edge from previous step
      // (it is actually top edge of some top
      // triangle in row j - 1)
      triang_edges_.row(it) << curr_e_ - Ni_ * 7 + i + 2, curr_e_ - 2, curr_e_ - 1;
      return;
    }
    // in general we do same thing as
    // in case (j == 1), but as we
    // have put more bottom edges in mesh
    // then in case (j == 1), we should
    // now decrement curr_e_ in other way
    triang_edges_.row(it) << curr_e_ - Ni_ * 6 + 1, curr_e_ - 2, curr_e_ - 1;
  }
  void triang_r(idx i, idx j) {
    idx it = (j * Ni_ + i) * 4 + 1;
    idx it0 = (i == Ni_ - 1) ? static_cast<idx>(boundaries::SOLID_WALL) :
      4 * (j * Ni_ + i + 1) + 3;
    idx it1 = (j * Ni_ + i) * 4 + 2;
    idx it2 = (j * Ni_ + i) * 4 + 0;
    triang_triangs_.row(it) << it0, it1, it2;
    idx ip0 = (Ni_ + 1) * j + i + 1;
    idx ip1 = (Ni_ + 1) * (j + 1) + i + 1;
    idx ip2 = num_verts_ + Ni_ * j  + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    edge_points_.row(curr_e_) << ip0, ip1;
    edge_triangs_.row(curr_e_++) << it, it0;
    edge_points_.row(curr_e_) << ip1, ip2;
    edge_triangs_.row(curr_e_++) << it, it1;
    triang_edges_.row(it) << curr_e_ - 2, curr_e_ - 1, curr_e_ - 4;
  }
  void triang_t(idx i, idx j) {
    idx it = (j * Ni_ + i) * 4 + 2;
    idx it0 = (j == Nj_ - 1) ? static_cast<idx>(boundaries::SOLID_WALL) :
      4 * ((j + 1) * Ni_ + i) + 0;
    idx it1 = (j * Ni_ + i) * 4 + 3;
    idx it2 = (j * Ni_ + i) * 4 + 1;
    triang_triangs_.row(it) << it0, it1, it2;
    idx ip0 = (Ni_ + 1) * (j + 1) + i + 1;
    idx ip1 = (Ni_ + 1) * (j + 1) + i;
    idx ip2 = num_verts_ + Ni_ * j + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    edge_points_.row(curr_e_) << ip1, ip0;
    edge_triangs_.row(curr_e_++) << it, it0;
    edge_points_.row(curr_e_) << ip1, ip2;
    edge_triangs_.row(curr_e_++) << it, it1;
    triang_edges_.row(it) << curr_e_ - 2, curr_e_ - 1, curr_e_ - 3;
  }
  void triang_l(idx i, idx j) {
    idx it = (j * Ni_ + i) * 4 + 3;
    idx it0 = (i == 0) ? static_cast<idx>(boundaries::SOLID_WALL) :
      4 * (j * Ni_ + i - 1) + 1;
    idx it1 = (j * Ni_ + i) * 4 + 0;
    idx it2 = (j * Ni_ + i) * 4 + 2;
    triang_triangs_.row(it) << it0, it1, it2;
    idx ip0 = (Ni_ + 1) * (j + 1) + i;
    idx ip1 = (Ni_ + 1) * j + i;
    idx ip2 = num_verts_ + Ni_ * j + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    if (i == 0) {
      // we add left edge of this square-box
      edge_points_.row(curr_e_) << ip1, ip0;
      edge_triangs_.row(curr_e_++) << it, it0;
      // so now our left edge is exactly
      // curr_e_-1
      triang_edges_.row(it) << curr_e_ - 1, curr_e_ - 6, curr_e_ - 2;
      return;
    }
    if (i == 1) {
      // but in this case we already have edge
      // with this coordinates
      // it's right edge of square {0, j}
      // so we should decrement curr_e_ in other way
      if (j == 0)
        triang_edges_.row(it) << curr_e_ - 12, curr_e_ - 5, curr_e_ - 1;
      else
        triang_edges_.row(it) << curr_e_ - 11, curr_e_ - 5, curr_e_ - 1;
      return;
    }
    // in general we do same thing as
    // in case (i == 1), but as the 
    // previous square had number {i-1, j}
    // and i-1 != 0, we now should decrement
    // curr_e_ on smaller value then in case
    // (i == 1)
    if (j == 0)
      triang_edges_.row(it) << curr_e_ - 11, curr_e_ - 5, curr_e_ - 1;
    else
      triang_edges_.row(it) << curr_e_ - 10, curr_e_ - 5, curr_e_ - 1;
  }
};

void testAll();

int main() {
  try {
    testAll();
    return 0;
  }
  catch (const ParserError& e) {
    std::cerr << "Parser error:" << std::endl;
    throw e;
  }
  catch (const DimensionError& e) {
    std::cerr << "Dimension manager error:" << std::endl;
    throw e;
  }
  catch (const MeshError& e) {
    std::cerr << "Mesh error:" << std::endl;
    throw e;
  }
  catch (const SolverError& e) {
    std::cerr << "Solver error:" << std::endl;
    throw e;
  }
  catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unrecognized exception type " << std::endl;
  }
  return -1;
}

void testValueFields() {
  using Eigen::Array2d;
  using Eigen::Array3d;

  StructTriangMesh M{ 2, 2, 0.1 };
  Bathymetry B{ std::move(M) };

	std::cout << B.str_size() << '\n';
  for (size_t i = 0; i < B.str_size(); i++)
    B.at_node(i) = -1.;

  VolumeField Uv{ B, B.mesh().num_triangles() };
	std::cout << Uv.storage_size() << '\n';
  Uv.prim(0) = Array<3>{ 0., 0., 15. };
  Uv.cons(1) = Array<3>{ 2., 0., 15. };
  std::cout << Uv.prim(1).get() << std::endl;

  EdgeField Ue{ B, B.mesh().num_triangles() * 3 };
  Ue.cons(0, 0, 1) = Array<3>{ 3., 0., 15. };
  Ue.vel (0, 0, 1) = Array<2>{ 2., 3. };
  std::cout << Ue.prim(0, 0, 1).get() << std::endl;
}

void testGaussWave() {

  using Eigen::Array3d;
  StructTriangMesh m{ 80, 80, 0.1 };
  Bathymetry b{ std::move(m) };

  for (size_t i = 0; i < b.str_size(); i++) {
    b.at_node(i) = 0.;
  }

  VolumeField v0{ b, b.mesh().num_triangles() };
  for (size_t i = 0; i < b.mesh().num_triangles(); ++i) {
    double r = (b.mesh().t(i) - Point{4., 4.}).matrix().squaredNorm();
    v0.prim(i) = Array3d{ 1. + exp(-5. * r), 0., 0. };
  }

	KurganovSpaceDisc disc{std::move(b), std::move(v0)};
	TimeDisc<KurganovSpaceDisc> solv{ &disc, EulerSolver<KurganovSpaceDisc>};

	double dt = 0.001;
	disc.Dump("out0.dat");
	solv.Step(dt);
	disc.Dump("out1.dat");
}

void testSimpleRunup() {

	using Eigen::Array3d;

  StructTriangMesh mesh{ 40, 40, 0.1 };
  Bathymetry bathymetry{ std::move(mesh) };

	Parser parser("config.ini");
	DimensionManager dimer(parser);
	BallTest test(parser, dimer, bathymetry.mesh());

  for (size_t i = 0; i < bathymetry.mesh().num_nodes(); i++) {
		const auto& p = bathymetry.mesh().p(i);
    bathymetry.at_node(i) = test.B(p[0], p[1]);
  }

  VolumeField v0{ bathymetry, bathymetry.mesh().num_triangles() };

	size_t i;
	#pragma omp parallel for
  for (i = 0; i < bathymetry.mesh().num_triangles(); ++i) {
 
		const auto& ip = bathymetry.mesh().triang_points(i);
		Point	p0 = bathymetry.mesh().p(ip[0]);
		Point	p1 = bathymetry.mesh().p(ip[1]);
		Point	p2 = bathymetry.mesh().p(ip[2]);

		double hc = triang_average(p0, p1, p2, [&](const Point& p) { return test.h(p[0], p[1], 0.); } );
		double uc = triang_average(p0, p1, p2, [&](const Point& p) { return test.u(p[0], p[1], 0.); } );
		double vc = triang_average(p0, p1, p2, [&](const Point& p) { return test.v(p[0], p[1], 0.); } );

		v0.prim(i) = Array3d{ hc + bathymetry.at_cell(i), uc, vc};
  }

	KurganovSpaceDisc disc{std::move(bathymetry), std::move(v0)};

	disc.Dump("out_0.dat");
	disc.reconstruct_all();
	disc.compute_fluxes();
	disc.dump("f_0.dat");

	TimeDisc<KurganovSpaceDisc> solv{ &disc, EulerSolver<KurganovSpaceDisc>};
	double t_end = 2.;

	int limiter = 0, max_nb_steps = 500;

	for (double t = 0., dt = 1e-3; t <= t_end; t += (dt = solv.CFLdt())) {
		printf("\rStep #%d:  t=%.5e,  dt = %.5e", ++limiter, t, dt);
		solv.Step(dt);
		if (limiter > max_nb_steps) {
			break;
		}
	}
	
	printf("\n");	
	disc.Dump("out_2.dat");
	disc.reconstruct_all();
	disc.compute_fluxes();
	disc.dump("f_2.dat");
}

void testAll() {
	testValueFields();
  //testGaussWave();
	testSimpleRunup();
}
