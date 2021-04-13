#pragma once
#include <Solver.h>

struct KurganovSpaceDisc : BaseSpaceDisc<KurganovSpaceDisc> {
	using BaseSpaceDisc::BaseSpaceDisc;
	Array<3> UserFlux(idx e, idx from, idx to);
};
