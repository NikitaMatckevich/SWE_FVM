#pragma once
#include <Solver.h>
 
struct EulerTimeSolver : BaseTimeSolver {

	using BaseTimeSolver::BaseTimeSolver;

	virtual void UserStep(double dt) override {
  	const auto& m = SpDisc_->mesh();
  	for (size_t i = 0; i < m.num_triangles(); ++i)
    	SpDisc_->getVolField().cons(i) += dt * SpDisc_->RHS(i);
	}
};


