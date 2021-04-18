#pragma once
#include "Solver.h"

struct KurganovSpaceDisc : BaseSpaceDisc<KurganovSpaceDisc> {
	using BaseSpaceDisc::BaseSpaceDisc;
	Array<3> UserFlux(Idx e, Idx from, Idx to);
};

