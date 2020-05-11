#include "Solver.h"

namespace utils {

  BaseFluxMethod::BaseFluxMethod(VolumeField&& v0) : MUSCLObject(std::move(v0)) {}

  using Array = BaseFluxMethod::Array;

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
} // namespace utils
