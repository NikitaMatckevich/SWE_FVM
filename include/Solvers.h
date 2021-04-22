#pragma once
#include <TimeDisc.h>

void EulerSolver(TimeDisc * const td, double dt);

void SSPRK3Solver(TimeDisc * const td, double dt);
