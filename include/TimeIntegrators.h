#pragma once
#include <Solver.h>

template <class SpaceDisc>
void EulerSolver(TimeDisc<SpaceDisc> * const td, double dt) {
	SpaceDisc* sd = td->getSpaceDisc();
	const auto& m = sd->mesh();
	for (size_t i = 0; i < m.num_triangles(); ++i) {
		sd->getVolField().cons(i) += td->RHS(i, dt);
	}
}


