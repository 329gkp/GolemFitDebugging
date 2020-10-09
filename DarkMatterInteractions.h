#ifndef _H_DMnu_
#define _H_DMnu_

#include <iostream>
#include <cmath>
#include <SQuIDS/const.h>
#include <nuSQuIDS/tools.h>

namespace DMnu {

  enum class Interaction {SS, SF, FV, FS};
  // SS : scalar dark matter, scalar mediator
  // SF : scalar dark matter, fermion mediator
  // FV : fermion dark matter, vector mediator
  // FS : fermion dark matter, scalar mediator

  struct ModelParameters {
    const Interaction type;
    const double g;
    const double mx;
    const double mphi;
    ModelParameters(Interaction type, double g, double mx, double mphi):
      type(type),g(g),mx(mx),mphi(mphi){};
  };

  double TotalCrossSection(double enu, ModelParameters mp);
  double DifferenrtialCrossSection(double enu_i, double enu_f, ModelParameters mp);

  double DarkMatterDensity(double b, double l, double x);
  double DarkMatterColumnDensity(double b, double l);

} // close namespace
#endif
