#include <Solvers.h>

void EulerSolver(TimeDisc * const td, double dt) { 
	SpaceDisc* sd = td->GetSpaceDisc();
	const auto& m = sd->Mesh();

	for (size_t i = 0; i < m.NumTriangles(); ++i) {
		sd->GetVolField().cons(i) += td->RHS(i, dt);
	}
}
