#pragma once
#include <TimeDisc.h>

namespace Solvers {

void Euler(TimeDisc * const td, double dt);
void SSPRK3(TimeDisc * const td, double dt);

}
