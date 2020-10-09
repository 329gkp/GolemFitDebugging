#ifndef _EARTH_H_INCLUDED
#define _EARTH_H_INCLUDED

#include <iostream>
#include <cmath>
#include <nuSQuIDS/tools.h>

namespace golemfit{
  namespace Earth {
    double RhoEarth(double theta, double x);
    double ColumnDensity(double theta);
    const double rE = 6371.; // km
    const double d = 1.948; // center of IC
  } // close namespace
} // close namespace
#endif
