#include "DarkMatterInteractions.h"

#define SQ(x)      ((x)*(x))                        // x^2
#define CUBE(x)      ((x)*(x)*(x))                        // x^3

namespace {

double costh(double enu_i, double enu_f, double mx){
  return mx/(1./enu_f - 1./enu_i);
}

double dcosthdE(double enu_f, double mx){
  return mx/SQ(enu_f);
}

double dsigmaSS(double enu_i, double enu_f, double g, double mx, double mphi){
  double x = costh(enu_i,enu_f,mx);
  double num = (1.-x)*SQ(enu_i)*mx;
  double denom = ((1.-x)*enu_i+mx)*SQ((1-x)*enu_i*SQ(mphi)+mx*(SQ(mphi)-2.*(x-1)*SQ(enu_i)));
  return dcosthdE(enu_f,mx)*(SQ(g)/(16.*M_PI))*(num/denom);
}

double sigmaSS(double enu, double g, double mx, double mphi){
  double C = 4.*SQ(enu)*mx+2.*enu*SQ(mphi)+mx*SQ(mphi);
  double B = C*log((SQ(mphi)*(2.*enu+mx))/C);
  double A = 4.*SQ(enu)*mx;
  double num = A + B;
  double denom = SQ(enu)*SQ(mx)*C;
  double coef = -SQ(g)/(64.*M_PI);
  return coef*num/denom;
}

double dsigmaFS(double enu_i, double enu_f, double g, double mx, double mphi){
  double x = costh(enu_i,enu_f,mx);
  double num = (x-1)*SQ(enu_i*mx)*(2*(x-1)*enu_i*mx+(x-1*SQ(enu_i)-2.*SQ(mx)));
  double denom = SQ(mx-(x-1)*enu_i)*SQ((x-1)*enu_i*SQ(mphi)-mx*(SQ(mphi)-2.*(x-1)*SQ(enu_i)));
  return dcosthdE(enu_f,mx)*(SQ(g)/(8.*M_PI))*(num/denom);
}

double sigmaFS(double enu, double g, double mx, double mphi){
  double X = 4.*SQ(enu)*mx+2.*enu*SQ(mphi)+mx*SQ(mphi);

  double A = enu*mx;
  double B = enu*SQ(mx)/(2.*enu+mx);
  double C1 = enu*SQ(mx)*SQ(mphi)*(SQ(mphi)-4.*SQ(mx));
  double C2 = (2.*enu*mx+SQ(mphi))*X;
  double C = C1/C2;
  double D = (enu*mx*(SQ(mphi)-4.*SQ(mx)))/(2.*enu*mx+SQ(mphi));
  double E = (SQ(mphi)-2.*SQ(mx))*log((SQ(mphi)*(2.*enu+mx)/X));
  double F = -SQ(g)/(32.*M_PI*SQ(enu*mx));

  return F*(A-B-C+D+E);
}

double dsigmaFV(double enu_i, double enu_f, double g, double mx, double mphi){
  double x = costh(enu_i,enu_f,mx);
  double num = SQ(enu_i*mx)*((1.-SQ(x))*enu_i*mx + SQ(1-x)*SQ(enu_i) + (1+x)*SQ(mx));
  double denom = ((1-x)*enu_i+mx)*((1-x)*enu_i*SQ(mx) + mx*SQ(SQ(mphi)-2.*(x-1)*SQ(enu_i)));
  return dcosthdE(enu_f,mx)*(SQ(g)/(4.*M_PI))*(num/denom);
}

double sigmaFV(double enu, double g, double mx, double mphi){
  double X = mx*SQ(mphi) + 4.*SQ(enu)*mx + 2.*enu*SQ(mphi);
  double A = SQ(mphi) + SQ(mx) - 2*enu*mx;
  double B = SQ(mphi)*(2.*enu+mx);
  double C = 2.*enu*(4.*SQ(enu)*mx+enu*(SQ(mx)+2.*SQ(mphi))+mx*SQ(mphi));
  double D = (2.*enu+mx)*X;

  double coef=SQ(g)/(16.*M_PI*SQ(enu*mx));

  return -coef*(A*log(B/X) + 4.*SQ(enu)*(-1.-SQ(mx)/SQ(mphi)+C/D));
}

double dsigmaSF(double enu_i, double enu_f, double g, double mx, double mphi){
  double x = costh(enu_i,enu_f,mx);
  double num = SQ(g)*(1+x)*SQ(SQ(enu_i))*(mx*mx*mx*mx*mx)*SQ((1-x)*enu_i+2.*mx);
  double denom = SQ(SQ(mphi) - mx *(2*enu_i+mx))*CUBE((1-x)*enu_i + mx)*SQ(enu_i*((x-1)*SQ(mphi) - (x+1)*SQ(mx)) + CUBE(mx) - mx*SQ(mphi));
  return dcosthdE(enu_f,mx)*(num/denom/(4.*M_PI));
}

double sigmaSF(double enu, double g, double mx, double mphi){
  double A = (8.*SQ(enu)*mx)/((2.*enu+mx)*SQ(SQ(mphi)-mx*(2.*enu+mx)));
  double B = 4./(SQ(mx)-SQ(mphi)-2.*enu*mx);
  double C = 8./(mx*(2.*enu+mx)-SQ(mphi));
  //double D = (3./SQ(enu) + (4.*enu*mx-2.*(SQ(mx)+3.*SQ(mphi)))/(enu*mx*(mx*(2.*enu+mx)-SQ(mphi))))*log(1+(4.*SQ(enu)*mx)/(SQ(mx)*(2.*enu+mx)-CUBE(mx)));
  double D1 = 1+(4.*SQ(enu)*mx)/(SQ(mphi)*(2.*enu+mx)-CUBE(mx));
  double D2 = (3./SQ(enu) + (4.*enu*mx-2.*(SQ(mx)+3.*SQ(mphi)))/(enu*mx*(mx*(2.*enu+mx)-SQ(mphi))));
  double D = D2*log(D1);
  return (SQ(g)/(64.*M_PI))*(A+B+C+D);
}

} // unname namespace

namespace DMnu {
squids::Const units;

double TotalCrossSection(double enu, ModelParameters mp) {
  switch (mp.type) {
    case Interaction::SS:
      return sigmaSS(enu, mp.g, mp.mx, mp.mphi);
    case Interaction::SF:
      return sigmaSF(enu, mp.g, mp.mx, mp.mphi);
    case Interaction::FV:
      return sigmaFV(enu, mp.g, mp.mx, mp.mphi);
    case Interaction::FS:
      return sigmaFS(enu, mp.g, mp.mx, mp.mphi);
    default:
      return std::numeric_limits<double>::signaling_NaN();
  }
}

double DifferentialCrossSection(double enu_i, double enu_f, ModelParameters mp) {
  switch (mp.type) {
    case Interaction::SS:
      return dsigmaSS(enu_i, enu_f, mp.g, mp.mx, mp.mphi);
    case Interaction::SF:
      return dsigmaSF(enu_i, enu_f, mp.g, mp.mx, mp.mphi);
    case Interaction::FV:
      return dsigmaFV(enu_i, enu_f, mp.g, mp.mx, mp.mphi);
    case Interaction::FS:
      return dsigmaFS(enu_i, enu_f, mp.g, mp.mx, mp.mphi);
    default:
      return std::numeric_limits<double>::signaling_NaN();
  }
}

double DarkMatterDensity(double b, double l, double x){
  double R = 8.5; // solar radius,  kpc
  double r = sqrt(x*x+R*R-2.*R*x*cos(b)*cos(l));
  // einasto profile
  double a = 0.17;
  double rs = 25.7;
  double rhos = 0.4 * exp((2./a)*(pow(R/rs,a)-1.));
  return rhos*exp((-2/a)*(pow(r/rs,a)-1.));
}

double DarkMatterColumnDensity(double b, double l) {
  nusquids::IntegrateWorkspace ws = nusquids::IntegrateWorkspace(5000);
  double xfinal = 40.; // kpc
  double u = units.GeV*(1.e3*units.parsec)/(units.cm*units.cm*units.cm);
  return integrate(ws,[&](double x) { return DarkMatterDensity(b,l,xfinal-x)*u;},0,xfinal, 1.e-5);
}

} // close namespace
