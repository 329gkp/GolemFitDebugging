/*
 * GolemParameters.cpp
 *
 *  Created on: Nov 21, 2017
 *      Author: aschneider
 */

#include <LeptonWeighter/ParticleType.h>
#include "GolemParameters.h"
#include "Event.h"

namespace golemfit {

/*************************************************************************************************************
 * METAGOLEM
 * **********************************************************************************************************/

std::ostream& operator<<(std::ostream& os, const FitParameters& fp){
  os << "Conventional normalization " << fp.convNorm << std::endl;;
  os << "Prompt normalization " << fp.promptNorm << std::endl;
  os << "Atmospheric muon normalization " << fp.muonNorm<< std::endl;
  os << "Astro component overall normalization " << fp.astroNorm<< std::endl;

  os << "Astro E component " << fp.getEFlavorRatio()<< std::endl;
  os << "Astro Mu component " << fp.getMuFlavorRatio()<< std::endl;
  os << "Astro Tau component " << fp.getTauFlavorRatio()<< std::endl;

  os << "Astroparticle balance " << fp.astroParticleBalance<< std::endl;
  os << "Astro gamma " << fp.astroDeltaGamma<< std::endl;
  os << "Astro cutoff " << fp.cutoffEnergy<< std::endl;
  os << "CR gamma " << fp.CRDeltaGamma<< std::endl;
  os << "Conv pi/k ratio " << fp.piKRatio<< std::endl;
  os << "Conv particle balance " << fp.NeutrinoAntineutrinoRatio<< std::endl;
  os << "Conv zenith correction " << fp.zenithCorrection<< std::endl;
  os << "Dark normalization " << fp.darkNorm<< std::endl;
  os << "DOM efficiency " << fp.domEfficiency<< std::endl;
  os << "Astro component overall normalization second " << fp.astroNormSec << std::endl;
  os << "Astro gamma second " << fp.astroDeltaGammaSec << std::endl;
  os << "Holeice forward " << fp.holeiceForward << std::endl;
  os << "Anisotropy scale " << fp.anisotropyScale << std::endl;
  os << "WP " << fp.barrWP << std::endl;
  os << "WM " << fp.barrWM << std::endl;
  os << "ZP " << fp.barrZP << std::endl;
  os << "ZM " << fp.barrZM << std::endl;
  os << "YP " << fp.barrYP << std::endl;
  os << "YM " << fp.barrYM << std::endl;
  os << "HP " << fp.barrHP << std::endl;
  os << "HM " << fp.barrHM << std::endl;
  return os;
}

unsigned int null_categorizer(const Event& e) {
  return 0;
}

unsigned int symmetric_categorizer(const Event& e) {
  return 0;
}

bool check_fit_parameters_equality(std::vector<FitParameters> fp0, FitParameters fp1, FitParametersFlag fpf) {
  bool result = true;
  for(auto fp : fp0)
    result *= check_fit_parameters_equality(fp,fp1,fpf);
  return result;
}
bool check_fit_parameters_equality(FitParameters fp0, FitParameters fp1, FitParametersFlag fpf) {
  return ((fpf.convNorm && (fp0.convNorm == fp1.convNorm)) || !fpf.convNorm) && ((fpf.promptNorm && (fp0.promptNorm == fp1.promptNorm)) || !fpf.promptNorm) && ((fpf.muonNorm && (fp0.muonNorm == fp1.muonNorm)) || !fpf.muonNorm) && ((fpf.astroNorm && (fp0.astroNorm == fp1.astroNorm)) || !fpf.astroNorm) &&
      ((fpf.astroFlavorAngle1 && (fp0.astroFlavorAngle1 == fp1.astroFlavorAngle1)) || !fpf.astroFlavorAngle1) &&
      ((fpf.astroFlavorAngle2 && (fp0.astroFlavorAngle2 == fp1.astroFlavorAngle2)) || !fpf.astroFlavorAngle2) &&
      //((fpf.astroENorm && (fp0.astroENorm == fp1.astroENorm)) || !fpf.astroENorm) &&
      //((fpf.astroMuNorm && (fp0.astroMuNorm == fp1.astroMuNorm)) || !fpf.astroMuNorm) &&
      //((fpf.astroTauNorm && (fp0.astroTauNorm == fp1.astroTauNorm)) || !fpf.astroTauNorm) &&
      ((fpf.astroParticleBalance && (fp0.astroParticleBalance == fp1.astroParticleBalance)) || !fpf.astroParticleBalance) && ((fpf.astroDeltaGamma && (fp0.astroDeltaGamma == fp1.astroDeltaGamma)) || !fpf.astroDeltaGamma) && ((fpf.cutoffEnergy && (fp0.cutoffEnergy == fp1.cutoffEnergy)) || !fpf.cutoffEnergy) && ((fpf.CRDeltaGamma && (fp0.CRDeltaGamma == fp1.CRDeltaGamma)) || !fpf.CRDeltaGamma) && ((fpf.piKRatio && (fp0.piKRatio == fp1.piKRatio)) || !fpf.piKRatio) && ((fpf.NeutrinoAntineutrinoRatio && (fp0.NeutrinoAntineutrinoRatio == fp1.NeutrinoAntineutrinoRatio)) || !fpf.NeutrinoAntineutrinoRatio) && ((fpf.zenithCorrection && (fp0.zenithCorrection == fp1.zenithCorrection)) || !fpf.zenithCorrection) && ((fpf.darkNorm && (fp0.darkNorm == fp1.darkNorm)) || !fpf.darkNorm) && ((fpf.domEfficiency && (fp0.domEfficiency == fp1.domEfficiency)) || !fpf.domEfficiency) && ((fpf.holeiceForward && (fp0.holeiceForward == fp1.holeiceForward)) || !fpf.holeiceForward) && ((fpf.anisotropyScale && (fp0.anisotropyScale == fp1.anisotropyScale)) || !fpf.anisotropyScale) && ((fpf.astroNormSec && (fp0.astroNormSec == fp1.astroNormSec)) || !fpf.astroNormSec) && ((fpf.astroDeltaGammaSec && (fp0.astroDeltaGammaSec == fp1.astroDeltaGammaSec)) || !fpf.astroDeltaGammaSec);
}


} // namespace golemfit
