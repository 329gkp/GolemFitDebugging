#ifndef _GOLEMTOOLS_H
#define _GOLEMTOOLS_H

#include <functional>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>
#include <memory>
#include <fstream>
#include <LeptonWeighter/ParticleType.h>
#include <nuSQuIDS/xsections.h>

namespace golemfit {
  namespace tools {
    using namespace nusquids;

    bool isNeutrino(LW::ParticleType p);
    bool isAntineutrino(LW::ParticleType p);

    std::string CheckedFilePath(std::string FilePath, bool quiet = false);

    class ScaledNeutrinoCrossSections : public NeutrinoDISCrossSectionsFromTables {
    private:
      // Energy bin edges of scale factors
      std::vector<double> edges;
      // Scale factors of those bins
      std::vector<double> scales;

    public:
      virtual ~ScaledNeutrinoCrossSections(){}
      /// \brief Default construct with built-in tables
      ScaledNeutrinoCrossSections()
        : edges(0), scales(0) {}
      ScaledNeutrinoCrossSections(std::string path)
        : NeutrinoDISCrossSectionsFromTables(path), edges(0), scales(0) {}
      ScaledNeutrinoCrossSections(std::vector<double> e, std::vector<double> s)
        : edges(e), scales(s) {}
      ScaledNeutrinoCrossSections(std::string path, std::vector<double> e, std::vector<double> s)
        : NeutrinoDISCrossSectionsFromTables(path), edges(e), scales(s) {}

      double GetScale(double Enu) const {
        unsigned int nedges = edges.size();
        if (nedges == 0) return 1.;
        if (Enu < edges[0]) return 1.;
        if (Enu > edges[nedges-1]) return 1.;
        for (unsigned int i=1;i<nedges;i++) {
          if (Enu < edges[i]) {
            return scales[i-1];
          }
        }
        return std::numeric_limits<double>::quiet_NaN();
      }
      /// \brief Returns the total neutrino cross section
      /// \details Used to interpolate the total cross sections.
      double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override {
        return GetScale(Enu)*NeutrinoDISCrossSectionsFromTables::TotalCrossSection(Enu,flavor,neutype,current);
      }
      /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
      /// \param E1 Incident lepton energy.
      /// \param E2 Outgoing lepton energy.
      /// \param flavor Flavor index.
      /// \param neutype Can be either neutrino or antineutrino.
      /// \param current Can be either CC or NC.
      /// \return The cross section in cm^2 GeV^-1.
      double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override {
        return GetScale(E1)*NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection(E1,E2,flavor,neutype,current);
      }
    };

  }// closing tools namespace
}// close golemfit namespace

#endif
