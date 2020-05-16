#include "Solver.h"

using Array = MUSCLObject::Array;

Array ElemFlux(const Eigen::Vector2d& n, const Array& Ucons) {
  Array res(Array::Zero());
  double h = Ucons[0];
  if (is_wet(h)) {
    auto hvel = Ucons.tail<2>();
    double hveln = hvel.matrix().dot(n);
    res[0] = hveln;
    res.tail<2>() = (hveln / h) * hvel + (0.5 * h * h) * n.array();
  }
  return res;
}
