#include "Earth.h"

namespace golemfit{  
  namespace Earth {  
    double RhoEarth(double theta, double x) {
      /* theta: zenith in radians
         x: path along trajectory in km

         returns density in g/cm^3
      */
      
      // theta = angle from down vector (0 = direction north pole...if you're at IceCube)
      // piecewise polynomial fit to Reference earth model STW105
      // you could also load a Ref earth model if you want.
    
      double r = sqrt(pow(rE-d,2) + pow(x,2) + 2. * (rE-d) * x * cos(theta));
      double p1, p2, p3 = 0;
      if (r < 1221.) {
        p1 = -0.0002177;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
      }
      else if (r < 3480.) {
        p1 = -0.0002409;
        p2 = 0.1416;
        p3 = 1.234e+04;
      }
      else if (r < 5721.) {
        p1 = -3.764e-05;
        p2 = -0.1876;
        p3 = 6664.;
      }
      else if (r < 5961.) {
        p1 = 0.;
        p2 = -1.269;
        p3 = 1.131e+04;
      }
      else if (r < 6347.) {
        p1 = 0.;
        p2 = -.725;
        p3 = 7887.;
      }
      else if (r < 6356.) {
        p1 = 0.;
        p2 = 0.;
        p3 = 2900.;
      }
      else if (r < 6368.) {
        p1 = 0.;
        p2 = 0.;
        p3 = 2600.;
      }
      else {
        p1 = 0.;
        p2 = 0.;
        p3 = 1020.;
      }

      double rho = p1 * r*r + p2 * r + p3;

      return rho*1.0e-3; // g/cm^3.1.0e-3 conversion factor from kg/m^3 to g/cm^3
    }
  
    double ColumnDensity(double theta) {
      /* theta: zenith in radians

         returns integrated density in g/cm^2
      */
      double xmax = sqrt(pow(rE-d, 2.)*pow(cos(theta), 2.)+d*(2*rE-d))-(rE-d)*cos(theta);

      double kmTocm = 1.0e5;
      nusquids::IntegrateWorkspace ws = nusquids::IntegrateWorkspace(5000);
      return nusquids::integrate(ws,[&](double x) { return RhoEarth(theta,xmax-x);},0,xmax, 1.e-5) * kmTocm;
    }
  } // close namespace
} // close namespace
