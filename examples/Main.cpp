#include <DimensionManager.h>
#include <StructTriangMesh.h>
#include <Fluxes.h>
#include <Solvers.h>
#include <Tests.h>

#include <cmath>
#include <iostream>
#include <iomanip>

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
    //#ifdef __GNUC__
     // PrintStacktrace();
    //#endif
    throw e;
  }
  catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    throw e;
  }
  catch (...) {
    std::cerr << "Unrecognized exception type " << std::endl;
  }
  return -1;
}

/*
Array<3> OutputL2Errors(const TriangMesh& m, const Storage<3>& err, const std::string& filename) {

  std::ofstream fout(filename);
  fout << "h\thu\thv\n";

  Array<3> sum = Array<3>::Zero();

  for (size_t i = 0; i < m.NumTriangles(); ++i) {
    sum += m.Area(i) * err.col(i).square();
    fout << err(0, i) << '\t' << err(1, i) << '\t' << err(2, i) << '\n';
  }

  fout.close();
  return sum.sqrt();
}
*/

void dumpFields(const SpaceDisc& sd, const std::string& filename) {
  std::ofstream fout(filename);
  const auto& vol = sd.GetVolField();
  const auto& m = sd.GetDomain().GetTopology();
  for (size_t i = 0; i < m.NumTriangles(); i++) {
    fout << vol.w(i) << '\t' << vol.hu(i) << '\t' << vol.hv(i) << '\n';
  }
  fout.close();
}

/*
void DumpIntegrals(double t, const Array<3>& ints, const std::string& filename) {
  std::ofstream fout(filename, std::ios_base::app);
  fout << t << ints.transpose() << ' ' << ints[1] + ints[2] << '\n';
  fout.close();
}

void DumpAnalytic(const Test& test, const TriangMesh& m, double t, const std::string& filename) {

  std::ofstream fout(filename);
	fout << "x\ty\th\tw\tu\tv\thu\thv\tb\n";

	for (size_t i = 0; i < m.NumTriangles(); i++) {
		
		const auto& p = m.T(i);
    const double x = p[0], y = p[1];	
	
		fout << x << '\t' << y << '\t';
		fout << test.h(x,y,t) << '\t' << test.w(x,y,t) << '\t';
		fout << test.u(x,y,t) << '\t' << test.v(x,y,t) << '\t'; 
		fout << test.hu(x,y,t) << '\t' << test.hv(x,y,t) << '\t';
		fout << test.b(x,y) << '\n';
	}

	fout.close();
}

void TestValueFields() {
 
  StructTriangMesh m{ 2, 2, 0.1 };
  Domain b{ std::move(m) };

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
*/

void dumpGeometry(const Domain& b, const std::string& geometryFile) {
  std::ofstream fout(geometryFile);
  const auto& m = b.GetTopology();
  for (size_t i = 0; i < m.NumTriangles(); i++) {
    const auto& t = b.T(i);
    fout << t[0] << '\t' << t[1] << '\t' << b.Z(i, t[0], t[1]) << '\n';
  }
  fout.close();  
}

void dumpGeometry2(const Domain& b, const std::string& geometryFile) {
  std::ofstream fout(geometryFile);
  const auto& m = b.Mesh();
  for (size_t i = 0; i < m.NumNodes(); i++) {
    const auto& p = m.P(i);
    fout << p[0] << '\t' << p[1] << '\t' << b.AtNode(i) << '\n';
  }
  fout.close();  
}

void dumpTopology(const Domain& b, const std::string& topologyFile) {
  std::ofstream fout(topologyFile);
  const auto& m = b.Mesh();
  fout << m.NumEdges() << '\n';
  for (size_t i = 0; i < m.NumEdges(); i++) {
    const auto& ep = m.EdgePoints(i);
    const auto& et = m.EdgeTriangs(i);
    fout << ep[0] << '\t' << ep[1] << '\t' << et[0] << '\t' << et[1] << '\n';
  }
  fout << m.NumTriangles() << '\n';
  for (size_t i = 0; i < m.NumTriangles(); i++) {
    const auto& tp = m.TriangPoints(i);
    const auto& te = m.TriangEdges(i);
    const auto& tt = m.TriangTriangs(i);
    fout << tp[0] << '\t' << tp[1] << '\t' << tp[2] << '\t';
    fout << te[0] << '\t' << te[1] << '\t' << te[2] << '\t';
    fout << tt[0] << '\t' << tt[1] << '\t' << tt[2] << '\n';
  }
  fout.close();
}

void dumpDomain(const Domain& b, const std::string& geometryFile, const std::string& topologyFile) {
  dumpGeometry2(b, geometryFile);
  dumpTopology(b, topologyFile);
}

void testGaussWave() {

  TriangMesh m("../examples/bowl.msh" );
  Domain b(&m);
  for (size_t i = 0; i < b.Size(); i++) {
    b.AtNode(i) = 0.;
  }

  dumpDomain(b, "geometry.dat", "topology.dat");

  VolumeField v0{ b, b.Mesh().NumTriangles() };
  for (size_t i = 0; i < b.Mesh().NumTriangles(); ++i) {
    double r = (b.Mesh().T(i) - Point{4., 4.}).matrix().squaredNorm();
    v0.prim(i) = Array<3>{ 1. + exp(-5. * r), 0., 0. };
  }

  SpaceDisc sd{Fluxes::HLL<Wavespeeds::Einfeldt>, std::move(b), std::move(v0)};
  TimeDisc td{&sd};

  double dt = 0.001;
  dumpFields(sd, "out0.dat");
  Solvers::Euler(&td, dt);
  dumpFields(sd, "out1.dat");
}

/*
Array<3> Process(const Test& test, size_t n, double l, double t_end, double cor = 0., double tau = 0.) {

  Domain bathymetry{StructTriangMesh{ n, n, l / n }};

  for (size_t i = 0; i < bathymetry.Mesh().NumNodes(); i++) {
		Point p = bathymetry.Mesh().P(i);
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

		Array<3> x = TriangAverage<3,100>(p0, p1, p2, [&](const Point& p) {
      return Array<3>{ test.h(p[0], p[1], 0.), test.u(p[0], p[1], 0.), test.v(p[0], p[1], 0.) };
    });
    x[0] += v0.b(i);
		v0.prim(i) = x;
  }

	SpaceDisc sd{Fluxes::HLLC<Wavespeeds::Einfeldt>, std::move(bathymetry), std::move(v0), cor, tau};
	TimeDisc td{&sd};

	int cnt = 0;
  int cnt_out = 0;
  double dt_out = 0.25 * t_end;
  double t_out = dt_out;

  Array<3> errL2 = Array<3>::Zero();

  std::string intfile = "int_" + std::to_string(n) + ".dat";

	for (double t = 0., dt = 1e-3; cnt_out < 4; t += dt) {

		printf("\rStep #%d:  t=%.5e,  dt = %.5e", ++cnt, t, dt);

		Solvers::SSPRK3(&td, dt);

    if (t >= t_out) {
      std::string outfile = "out_" + std::to_string(n) + "_" + std::to_string(cnt) + ".dat";
      std::string errfile = "err_" + std::to_string(n) + "_" + std::to_string(cnt) + ".dat";
      DumpFields(sd, outfile);
      errL2 = OutputL2Errors(sd.Mesh(), sd.CompareWith(test, t_out), errfile);
      t_out += dt_out;
      cnt_out++;
    }

    DumpIntegrals(t, sd.ComputeIntegrals(), intfile);
	}

  printf("\n");
  return errL2;
}

void TestConvergence() {

  const double l = 4.;
  
  Parser parser("config.ini");
	DimensionManager dimer(parser);
  ClassicThackerTest test(parser, dimer, 0.5*l, 0.5*l);
  const double cor = dimer.Unscale<Scales::source>(parser.Get("Common", "cor"));  
  const double tau = dimer.Unscale<Scales::source>(parser.Get("Common", "tau"));  
  const double t_end = M_PI / sqrt(8 + cor * cor);

  Array<3> err_p = Array<3>::Zero();

  std::ofstream fout("convergence_hllce.dat");
  fout << "n\th\thu\thv\tord_h\tord_hu\tord_hv\n";

  for (int k = 5; k <= 5; ++k) {
    size_t n = 1 << k;
    printf("Testing for n=%lu\n", n);
    Array<3> err_n = Process(test, n, l, t_end, cor);
    Array<3> ord_n = (err_p / err_n).log() / std::log(2);
    fout << 4 * n * n << ' ';
    fout << std::scientific << std::setprecision(6) << err_n[0] << ' ' << err_n[1] << ' ' << err_n[2] << ' ';
    fout << std::scientific << std::setprecision(6) << ord_n[0] << ' ' << ord_n[1] << ' ' << ord_n[2] << '\n';
    err_p = err_n;
  }

  fout.close();
}

void AddAnalyticPictures() {
  const double l = 4.;
  
  Parser parser("config.ini");
	DimensionManager dimer(parser);
  ClassicThackerTest test(parser, dimer, 0.5*l, 0.5*l);
  const double cor = dimer.Unscale<Scales::source>(parser.Get("Common", "cor"));  
  const double tau = dimer.Unscale<Scales::source>(parser.Get("Common", "tau"));  
  const double t_end = 2. * M_PI / sqrt(8 + cor * cor);

  size_t n = 128;
  StructTriangMesh m{ n, n, l / n };

  for (double t = 0; t <= t_end + 0.01; t += 0.25 * t_end) {
    std::string anafile = "ana_" + std::to_string(t) + ".dat";
    DumpAnalytic(test, m, t, anafile);
  }
}

void TestLakeAtRest() {

  const size_t n = 32;
  const double l = 4;
  Parser parser("config.ini");
	DimensionManager dimer(parser);
  LakeAtRestTest test(parser, dimer, 0.5*l, 0.5*l);
  Domain bathymetry{StructTriangMesh{ n, n, l / n }};

  for (size_t i = 0; i < bathymetry.Mesh().NumNodes(); i++) {
		Point p = bathymetry.Mesh().P(i);
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

		Array<3> x = TriangAverage<3,100>(p0, p1, p2, [&](const Point& p) {
      return Array<3>{ std::max(0., bathymetry.AtPoint(i,p)), test.u(p[0], p[1], 0.), test.v(p[0], p[1], 0.) };
    });
		v0.prim(i) = x;
  }

	SpaceDisc sd{Fluxes::HLLC<Wavespeeds::Einfeldt>, std::move(bathymetry), std::move(v0)};
	DumpFields(sd, "init.dat");
  TimeDisc td{&sd};

  double t_end = 0.5;
	int cnt = 0;
  int cnt_out = 0;
  double dt_out = 0.25 * t_end;
  double t_out = dt_out;

  Array<3> errL2 = Array<3>::Zero();

  std::string intfile = "Lint_" + std::to_string(n) + ".dat";

	for (double t = 0., dt = 1e-3; cnt_out < 4; t += dt) {

		printf("\rStep #%d:  t=%.5e,  dt = %.5e", ++cnt, t, dt);

		Solvers::SSPRK3(&td, dt);

    if (t >= t_out) {
      std::string outfile = "Lout_" + std::to_string(n) + "_" + std::to_string(cnt) + ".dat";
      std::string errfile = "Lerr_" + std::to_string(n) + "_" + std::to_string(cnt) + ".dat";
      DumpFields(sd, outfile);
      errL2 = OutputL2Errors(sd.Mesh(), sd.CompareWith(test, t_out), errfile);
      t_out += dt_out;
      cnt_out++;
    }

    DumpIntegrals(t, sd.ComputeIntegrals(), intfile);
	}

  printf("\n");
  std::cout << errL2 << std::endl; 
}
*/

void testAll() {
  testGaussWave();
	//testValueFields();
  //testGaussWave();
	//TestConvergence();
  //TestLakeAtRest();
  //AddAnalyticPictures();
}
