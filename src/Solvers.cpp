#include <Solvers.h>

void EulerSolver(TimeDisc * const td, double dt) { 
	SpaceDisc* sd = td->GetSpaceDisc();
	const auto& m = sd->Mesh();

  sd->ReconstructAll();
  sd->ComputeFluxes();
	for (size_t i = 0; i < m.NumTriangles(); ++i) {
		sd->GetVolField().cons(i) += td->RHS(i, dt);
	}
}

void SSPRK3Solver(TimeDisc * const td, double dt) {
	SpaceDisc* sd = td->GetSpaceDisc();
	const auto& m = sd->Mesh();

  sd->ReconstructAll();
  sd->ComputeFluxes();
  const VolumeField U0 = sd->GetVolField();
	for (size_t i = 0; i < m.NumTriangles(); ++i) {
		sd->GetVolField().cons(i) = U0.cons(i) + td->RHS(i, dt);
	}
  
  sd->ReconstructAll();
  sd->ComputeFluxes();
  const VolumeField U1 = sd->GetVolField();
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		sd->GetVolField().cons(i) = 0.75 * U0.cons(i) + 0.25 * U1.cons(i) + td->RHS(i, 0.25 * dt);
	}
  
  sd->ReconstructAll();
  sd->ComputeFluxes();
  const VolumeField& U2 = sd->GetVolField();
  for (size_t i = 0; i < m.NumTriangles(); ++i) {
		sd->GetVolField().cons(i) = (1./3.) * U0.cons(i) + (2./3.) * U2.cons(i) + td->RHS(i, (2./3.) * dt);
	}
}