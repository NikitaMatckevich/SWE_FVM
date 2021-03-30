#pragma once
#include <Solver.h>

struct KurganovSpaceDisc : BaseSpaceDisc {
	using BaseSpaceDisc::BaseSpaceDisc;
	Array UserFlux(idx e, idx from, idx to) override;
};
