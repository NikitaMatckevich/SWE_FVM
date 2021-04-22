#include <DimensionManager.h>
#include <StructTriangMesh.h>
#include <Fluxes.h>
#include <Solvers.h>
#include <Tests.h>

#include <omp.h>

#include <iostream>

void TestAll();

int main() {
  try {
    TestAll();
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

void TestValueFields() {
 
  StructTriangMesh m{ 2, 2, 0.1 };
  Bathymetry b{ std::move(m) };

	std::cout << b.Size() << '\n';
  for (size_t i = 0; i < b.Size(); i++)
    b.AtNode(i) = -1.;

  VolumeField uv{ b, b.Mesh().NumTriangles() };

	std::cout << uv.Size() << '\n';

  uv.prim(0) = Array<3>{ 0., 0., 15. };
  uv.cons(1) = Array<3>{ 2., 0., 15. };
  std::cout << uv.prim(1).Get() << std::endl;

  EdgeField ue{ b, b.Mesh().NumTriangles() * 3 };
  ue.cons(0, 0, 1) = Array<3>{ 3., 0., 15. };
  ue.vel (0, 0, 1) = Array<2>{ 2., 3. };
  std::cout << ue.prim(0, 0, 1).Get() << std::endl;
}

void TestGaussWave() {

  Bathymetry b{StructTriangMesh{ 80, 80, 0.1 }};

  for (size_t i = 0; i < b.Size(); i++) {
    b.AtNode(i) = 0.;
  }

  VolumeField v0{ b, b.Mesh().NumTriangles() };
  for (size_t i = 0; i < b.Mesh().NumTriangles(); ++i) {
    double r = (b.Mesh().T(i) - Point{4., 4.}).matrix().squaredNorm();
    v0.prim(i) = Array<3>{ 1. + exp(-5. * r), 0., 0. };
  }

	SpaceDisc sd{Fluxes::HLL<Wavespeeds::Einfeldt>, std::move(b), std::move(v0)};
	TimeDisc td{&sd};

	double dt = 0.001;
	sd.DumpFields("out0.dat");
	Solvers::Euler(&td, dt);
	sd.DumpFields("out1.dat");
}

void TestSimpleRunup() {

  Bathymetry bathymetry{StructTriangMesh{ 40, 40, 0.1 }};

	Parser parser("config.ini");
	DimensionManager dimer(parser);
	BallTest test(parser, dimer, bathymetry.Mesh());

  for (size_t i = 0; i < bathymetry.Mesh().NumNodes(); i++) {
		const auto& p = bathymetry.Mesh().P(i);
    bathymetry.AtNode(i) = test.b(p[0], p[1]);
  }

  VolumeField v0{ bathymetry, bathymetry.Mesh().NumTriangles() };

	size_t i;

	#pragma omp parallel for
  for (i = 0; i < bathymetry.Mesh().NumTriangles(); ++i) {
 
		const auto& ip = bathymetry.Mesh().TriangPoints(i);
		Point	p0 = bathymetry.Mesh().P(ip[0]);
		Point	p1 = bathymetry.Mesh().P(ip[1]);
		Point	p2 = bathymetry.Mesh().P(ip[2]);

		double hc = TriangAverage(p0, p1, p2, [&](const Point& p) { return test.h(p[0], p[1], 0.); } );
		double uc = TriangAverage(p0, p1, p2, [&](const Point& p) { return test.u(p[0], p[1], 0.); } );
		double vc = TriangAverage(p0, p1, p2, [&](const Point& p) { return test.v(p[0], p[1], 0.); } );

		v0.prim(i) = Array<3>{ hc + v0.b(i), uc, vc};
  }

	SpaceDisc sd{Fluxes::HLLC<Wavespeeds::Einfeldt>, std::move(bathymetry), std::move(v0)};
	sd.DumpFields("out_0.dat");

	TimeDisc td{&sd};
	double t_end = 2.;

	int limiter = 0, max_nb_steps = 10000;

	for (double t = 0., dt = 1e-3; t <= t_end; t += (dt = td.CFLdt())) {
		printf("\rStep #%d:  t=%.5e,  dt = %.5e", ++limiter, t, dt);
		Solvers::SSPRK3(&td, dt);
		if (limiter > max_nb_steps) {
			break;
		}
	}
	
	printf("\n");	
	sd.DumpFields("out_2.dat");
}

void TestAll() {
	//testValueFields();
  //testGaussWave();
	TestSimpleRunup();
}
