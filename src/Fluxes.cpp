#include <Fluxes.h>

namespace Wavespeeds {
 
  Array<2> Rusanov(double ul, double hl, double ur, double hr) {
    double cl = sqrt(hl);
    double cr = sqrt(hr);
    double aplus = std::max(abs(ul) + cl, abs(ur) + cr);
    return Array<2>::Constant(aplus);
  }

  Array<2> Davis(double ul, double hl, double ur, double hr) {
    double cl = sqrt(hl);
    double cr = sqrt(hr);
    return Array<2>{std::min(ul - cl, ur - cr), std::max(ul + cl, ur + cr)};
  }

  Array<2> Einfeldt(double ul, double hl, double ur, double hr) {
    double cl = sqrt(hl);
    double cr = sqrt(hr);
    
    double uRoe = (cl * ul + cl * ur) / (cl + cr);
    double cRoe = sqrt(0.5 * (hl + hr));

    return Array<2>{std::min(ul - cl, uRoe - cRoe), std::max(ur + cr, uRoe + cRoe)};
  }

}
