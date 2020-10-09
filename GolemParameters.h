/*
 * GolemParameters.h
 *
 *  Created on: Nov 21, 2017
 *      Author: aschneider
 */

#ifndef GOLEMPARAMETERS_H_
#define GOLEMPARAMETERS_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <boost/math/constants/constants.hpp>

#include "Event.h"
#include "GolemEnumDefinitions.h"
#include "DarkMatterInteractions.h"
#include "utils.h"

namespace golemfit{

struct FitParameters {
  float convNorm=1.0;
  float promptNorm=1.0;
  float CRDeltaGamma=0.0;
  float zenithCorrection=0.0;

  float muonNorm=1.0;
  float astroNorm=8.0/3.0;
  float astroFlavorAngle1=4./9.;
  float astroFlavorAngle2=0.0;

  //float astroENorm=1.0/3.;
  //float astroMuNorm=1.0/3.;
  //float astroTauNorm=1.0/3.;

  float astroParticleBalance=1.0;
  float astroDeltaGamma=0.0;
  float cutoffEnergy=10; // in log10
  float piKRatio=1.0;
  float NeutrinoAntineutrinoRatio=1.0;
  float darkNorm=1.e-20;
  float domEfficiency=0.99;
  float hqdomEfficiency=1.35;
  float holeiceForward=0;
  float anisotropyScale=1.0;
  float astroNormSec=1.0;
  float astroDeltaGammaSec=0.0;

  float barrWP = 0;
  float barrWM = 0;
  float barrYP = 0;
  float barrYM = 0;
  float barrZP = 0;
  float barrZM = 0;
  float barrHP = 0;
  float barrHM = 0;

  float icegrad0 = 0;
  float icegrad1 = 0;
  float icegrad2 = 0;

  float nuxs = 1.;
  float nubarxs = 1.;

  float kaonLosses = 0.0;

  FitParameters(){};
  FitParameters(sampleTag sample){
    if(sample == sampleTag::Sterile){
       astroDeltaGamma           = 0;
       astroDeltaGammaSec        = 0;
       // astroENorm                = 0;
       // astroMuNorm               = 0;
       astroNorm                 = 4.72/6.;
       astroParticleBalance      = 0;
       //astroTauNorm              = 0;
       cutoffEnergy              = 0;
       darkNorm                  = 0;
       muonNorm                  = 0;
       promptNorm                = 1.0;
       barrHM                    = 0;
       barrHP                    = 0;
       barrWM                    = 0;
       barrWP                    = 0;
       barrYM                    = 0;
       barrYP                    = 0;
       barrZM                    = 0;
       barrZP                    = 0;
       domEfficiency             = 1.27;
       holeiceForward            = -1.;
       NeutrinoAntineutrinoRatio = 1.;
       piKRatio                  = 1.;
       convNorm                  = 1.;
       CRDeltaGamma              = 0.;
       zenithCorrection          = 0.;
       hqdomEfficiency          = 1.35;
       nuxs = 1.;
       nubarxs = 1.;
       kaonLosses = 0.0;
    } else if (sample == sampleTag::HESE){
      darkNorm = 0;
      astroNormSec = 0;
    }
  }
  void setFlavorRatio(float f0, float f1, float f2) {
      float n = f0+f1+f2;
      f0 /= n;
      f1 /= n;
      f2 /= n;

      float a = golemfit::detail::convertFlavorRatioToAngle<0,float>()(f0, f1, f2);
      float b = golemfit::detail::convertFlavorRatioToAngle<1,float>()(f0, f1, f2);
      astroFlavorAngle1 = a;
      astroFlavorAngle2 = b;
  }
  float getEFlavorRatio() const {
    return (astroFlavorAngle2+1.0)*pow(astroFlavorAngle1, 0.5)/2.0;
  }
  float getMuFlavorRatio() const {
    return (-astroFlavorAngle2+1.0)*pow(astroFlavorAngle1, 0.5)/2.0;
  }
  float getTauFlavorRatio() const {
    return -pow(astroFlavorAngle1, 0.5) + 1.0;
  }
  FitParameters(sampleTag sample, FluxComponent component):FitParameters(sample){
    switch(component) {
      case FluxComponent::atmConv:
        convNorm=1.0;
        promptNorm=0.0;
        muonNorm=0.0;
        astroNorm=0.0;
        break;
      case FluxComponent::atmPrompt:
        convNorm=0.0;
        promptNorm=1.0;
        muonNorm=0.0;
        astroNorm=0.0;
        break;
      case FluxComponent::atmMuon:
        convNorm=0.0;
        promptNorm=0.0;
        muonNorm=1.0;
        astroNorm=0.0;
        break;
      case FluxComponent::diffuseAstro:
        convNorm=0.0;
        promptNorm=0.0;
        muonNorm=0.0;
        astroNorm=1.0;
        setFlavorRatio(1.0, 1.0, 1.0);
        break;
      case FluxComponent::diffuseAstro_e:
        convNorm=0.0;
        promptNorm=0.0;
        muonNorm=0.0;
        astroNorm=1.0;
        setFlavorRatio(1.0, 0.0, 0.0);
        break;
      case FluxComponent::diffuseAstro_mu:
        convNorm=0.0;
        promptNorm=0.0;
        muonNorm=0.0;
        astroNorm=1.0;
        setFlavorRatio(0.0, 1.0, 0.0);
        break;
      case FluxComponent::diffuseAstro_tau:
        convNorm=0.0;
        promptNorm=0.0;
        muonNorm=0.0;
        astroNorm=1.0;
        setFlavorRatio(0.0, 0.0, 1.0);
        break;
      default:
        throw std::runtime_error("Invalid flux component encountered");
    }
  }
};

std::ostream& operator<<(std::ostream& os, const FitParameters& fp);

// true means its fix in the minimization. CA.
// defaults corresponds to the ClassicHESE paramters
struct FitParametersFlag {
  bool convNorm=false;
  bool promptNorm=false;
  bool CRDeltaGamma=false;
  bool zenithCorrection=true;
  bool muonNorm=false;
  bool astroNorm=false;
  bool astroFlavorAngle1=true;
  bool astroFlavorAngle2=true;

  //bool astroENorm=true;
  //bool astroMuNorm=true;
  //bool astroTauNorm=true;

  bool astroParticleBalance=true;
  bool astroDeltaGamma=false;
  bool cutoffEnergy=true;
  bool piKRatio=false;
  bool NeutrinoAntineutrinoRatio=true;
  bool darkNorm=true;
  bool anisotropyScale=true;
  bool astroNormSec=true;
  bool astroDeltaGammaSec=true;

  bool barrWP=true;
  bool barrWM=true;
  bool barrZP=true;
  bool barrZM=true;
  bool barrYP=true;
  bool barrYM=true;
  bool barrHP=true;
  bool barrHM=true;

  bool domEfficiency=true;
  bool hqdomEfficiency=true;
  bool holeiceForward=true;

  bool icegrad0=true;
  bool icegrad1=true;
  bool icegrad2=true;

  bool nuxs=true;
  bool nubarxs=true;

  bool kaonLosses=true;

  FitParametersFlag(){};
  FitParametersFlag(bool flag):
    convNorm(flag),
    promptNorm(flag),
    CRDeltaGamma(flag),
    zenithCorrection(flag),

    muonNorm(flag),
    astroNorm(flag),
    astroFlavorAngle1(flag),
    astroFlavorAngle2(flag),

    //astroENorm(flag),
    //astroMuNorm(flag),
    //astroTauNorm(flag),

    astroParticleBalance(flag),
    astroDeltaGamma(flag),
    cutoffEnergy(flag),
    piKRatio(flag),
    NeutrinoAntineutrinoRatio(flag),
    darkNorm(flag),
    anisotropyScale(flag),
    astroNormSec(flag),
    astroDeltaGammaSec(flag),

    barrWP(flag),
    barrWM(flag),
    barrZP(flag),
    barrZM(flag),
    barrYP(flag),
    barrYM(flag),
    barrHP(flag),
    barrHM(flag),

    domEfficiency(flag),
    hqdomEfficiency(flag),
    holeiceForward(flag),

    icegrad0(flag),
    icegrad1(flag),
    icegrad2(flag),

    nuxs(flag),
    nubarxs(flag),

    kaonLosses(flag)
  {};
  FitParametersFlag(sampleTag sample):FitParametersFlag(true){
    if(sample == sampleTag::Sterile){
      // these are the default minimization parameters
      convNorm = false;
      promptNorm = false;
      CRDeltaGamma = false;
      astroNorm = false;
      astroDeltaGamma = false;
      //piKRatio = false;
      barrWP=false;
      barrWM=false;
      barrZP=false;
      barrZM=false;
      barrYP=false;
      barrYM=false;
      //barrHP=false;
      //barrHM=false;
      NeutrinoAntineutrinoRatio=false;
      domEfficiency=false;
      //hqdomEfficiency=true;
      holeiceForward=false;
      zenithCorrection=false;
      icegrad0 = false;
      icegrad1 = false;
      nuxs = false;
      nubarxs = false;
      kaonLosses = false;
    } else if (sample == sampleTag::HESE){
      convNorm = false;
      promptNorm = false;
      muonNorm = false;
      astroNorm = false;
      astroDeltaGamma = false;
      CRDeltaGamma = false;
      piKRatio = false;
    }
  }
};

bool check_fit_parameters_equality(FitParameters fp0, FitParameters fp1, FitParametersFlag fpf);
bool check_fit_parameters_equality(std::vector<FitParameters> fp0, FitParameters fp1, FitParametersFlag fpf);

struct Priors {
  float convNormCenter=1.;
  float convNormWidth=0.4;
  float muonNormCenter=1.;
  float muonNormWidth=0.5;
  float promptNormCenter=1.;
  float promptNormWidth=std::numeric_limits<float>::max();
  float crSlopeCenter=0.0; //
  float crSlopeWidth=0.01; //https://www.researchgate.net/profile/L_Marcelli/publication/265671393_Cosmic_Rays_under_the_Knee/links/5605240d08ae8e08c08ae3da/Cosmic-Rays-under-the-Knee.pdf
  float piKRatioCenter=1.0;
  float piKRatioWidth=0.1;
  float nuNubarRatioCenter=1.0;
  float nuNubarRatioWidth=0.1;
  float zenithCorrectionMultiplier=0.038;
  float domEffCenter = 0.99;
  float domEffWidth  = 100.0;
  float hqdomEffCenter = 1.35;
  float hqdomEffWidth  = 100.0;
  float holeiceForwardCenter = 0;
  float holeiceForwardWidth  = 0.5;

  float anisotropyScaleCenter=1.0;
  float anisotropyScaleWidth=0.2;
  float astroNormCenter= 8./3;
  float astroNormPower = 0.0;
  float astroNormWidth=std::numeric_limits<float>::infinity();
  float astroDeltaGammaCenter= 0.0;
  float astroDeltaGammaWidth=std::numeric_limits<float>::infinity();
  float astroFlavorAngle1Center=4./9.;
  float astroFlavorAngle1Width=std::numeric_limits<float>::infinity();
  float astroFlavorAngle2Center=0.;
  float astroFlavorAngle2Width=std::numeric_limits<float>::infinity();


  //float astroENormCenter=0.;
  //float astroENormWidth=std::numeric_limits<float>::max();
  //float astroMuNormCenter=0.;
  //float astroMuNormWidth=std::numeric_limits<float>::max();
  //float astroTauNormCenter=0.;
  //float astroTauNormWidth=std::numeric_limits<float>::max();

  float astroNormSecCenter=0.;
  float astroNormSecWidth=std::numeric_limits<float>::max();
  float astroDeltaGammaSecCenter=0.;
  float astroDeltaGammaSecWidth=std::numeric_limits<float>::max();
  //means n= 2.50, 4.72
  //cov=[[0.18527653, 0.94099016],
  //     [0.94099016, 7.59454058]]
  /*
  float astro1Comp2DNormCenter = 4.72/6.0;
  float astro1Comp2DNormWidth = sqrt(0.185276)/6.0;
  float astro1Comp2DDeltaGammaCenter = 0.0;
  float astro1Comp2DDeltaGammaWidth = sqrt(7.5945058);
  float astro1Comp2DCorrelation = 0.94099016/sqrt(0.185276*7.5945058)/6.;
  */

  float astro1Comp2DNormCenter           = 4.72/6;
  float astro1Comp2DNormWidth            = 0.36;
  float astro1Comp2DDeltaGammaCenter     = 0.00;
  float astro1Comp2DDeltaGammaWidth      = 0.36;
  float astro1Comp2DCorrelation          = 0.70;

  float barrWPCenter = 0;
  float barrWPWidth = 0.40;
  float barrWPMin = -0.5;
  float barrWPMax = 0.5;

  float barrWMCenter = 0;
  float barrWMWidth = 0.40;
  float barrWMMin = -0.5;
  float barrWMMax = 0.5;

  float barrYPCenter = 0;
  float barrYPWidth = 0.30;
  float barrYPMin = -0.5;
  float barrYPMax = 0.5;

  float barrYMCenter = 0;
  float barrYMWidth = 0.30;
  float barrYMMin = -0.5;
  float barrYMMax = 0.5;

  float barrZPCenter = 0;
  float barrZPWidth = 0.12;
  float barrZPMin = -0.2;
  float barrZPMax = 0.5;

  float barrZMCenter = 0;
  float barrZMWidth = 0.12;
  float barrZMMin = -0.25;
  float barrZMMax = 0.5;

  float barrHPCenter = 0;
  float barrHPWidth = 0.15;
  float barrHPMin = -0.5;
  float barrHPMax = 0.5;

  float barrHMCenter = 0;
  float barrHMWidth = 0.15;
  float barrHMMin = -0.5;
  float barrHMMax = 0.5;

  float icegrad0Center = 0;
  float icegrad0Width = 1.0;
  float icegrad0Max= 5;
  float icegrad0Min= -5;

  float icegrad1Center = 0;
  float icegrad1Width = 1.0;
  float icegrad1Max= 5;
  float icegrad1Min= -5;

  float icegrad2Center = 0;
  float icegrad2Width = 0.5;
  float icegrad2Max= 0.1;
  float icegrad2Min= -0.1;

  float icegrad01_correlation = 5.091035738186185100e-02;

  float zenithCorrectionCenter = 0.0;
  float zenithCorrectionWidth = 1.0;

  float nuXSCenter = 1.0;
  float nuXSWidth = 0.1;
  float nubarXSCenter = 1.0;
  float nubarXSWidth = 0.1;

  float kaonLossesCenter = 0.0;
  float kaonLossesWidth = 1.0;

  Priors(){}
  Priors(sampleTag sample){
    if(sample == sampleTag::Sterile){
      convNormCenter 	= 1;
      convNormWidth 	= 0.4;
      promptNormCenter 	= 1.;
      promptNormWidth 	= std::numeric_limits<float>::max();
      crSlopeCenter 	= 0.0;
      crSlopeWidth 	= 0.03; // This should go to 0.01
      piKRatioCenter 	= 1.0;  // will be remove for Barr
      piKRatioWidth 	= 0.1;  // will be remove for Barr
      nuNubarRatioCenter= 1.0;
      nuNubarRatioWidth = 0.025;
      domEffCenter 	= 1.27;
      domEffWidth 	= 0.123; // This is a 10% dom efficiency width. 31% of the DOM efficiency is the compensation factor.
      hqdomEffCenter 	= 1.35;
      hqdomEffWidth 	= 0.02; // This is a 10% dom efficiency width. 31% of the DOM efficiency is the compensation factor.
      holeiceForwardCenter = -1.0;
      holeiceForwardWidth  = 10.0; // We don't have a prior on this. We'll keep it flat.
      icegrad0Center = 0.0;
      icegrad0Width = 1.0;
      icegrad1Center = 0.0;
      icegrad1Width = 1.0;
      icegrad01_correlation = 5.091035738186185100e-02;
      nubarXSCenter     = 1.;
      nubarXSWidth      = 0.075;
      nubarXSCenter     = 1.;
      nubarXSWidth      = 0.03;
      kaonLossesCenter = 0.;
      kaonLossesWidth = 1.0;
    }
  }
  Priors(const FitParameters& fit_params){
   convNormCenter=fit_params.convNorm;
   muonNormCenter=fit_params.muonNorm;
   promptNormCenter=fit_params.promptNorm;
   crSlopeCenter= fit_params.CRDeltaGamma;
   piKRatioCenter=fit_params.piKRatio;
   nuNubarRatioCenter=fit_params.NeutrinoAntineutrinoRatio;
   domEffCenter=fit_params.domEfficiency;
   hqdomEffCenter=fit_params.hqdomEfficiency;
   holeiceForwardCenter=fit_params.holeiceForward;
   anisotropyScaleCenter=fit_params.anisotropyScale;
   astroNormCenter=fit_params.astroNorm;
   astroDeltaGammaCenter = fit_params.astroDeltaGamma;
   astroFlavorAngle1Center= fit_params.astroFlavorAngle1;
   astroFlavorAngle2Center= fit_params.astroFlavorAngle2;
   astroNormSecCenter = fit_params.astroNormSec;
   astroDeltaGammaSecCenter = fit_params.astroDeltaGammaSec;
   icegrad0Center = fit_params.icegrad0;
   icegrad1Center = fit_params.icegrad1;
   astro1Comp2DNormCenter = fit_params.astroNorm;
   astro1Comp2DDeltaGammaCenter = fit_params.astroDeltaGamma;
   nuXSCenter = fit_params.nuxs;
   nubarXSCenter = fit_params.nubarxs;
   kaonLossesCenter = fit_params.kaonLosses;
   //astro1Comp2DNormCenter = fit_params.astro1Comp2DNormWid;
   //astro1Comp2DDeltaGammaCenter = fit_params.astro1Comp2DDeltaGammaCenter;
  }
};

struct FitResult {
  FitParameters params;
  double likelihood;
  double aux_likelihood;
  unsigned int nEval, nGrad;
  bool succeeded;
  FitResult(){};
};

struct DataPaths {
  private:
   bool CheckDataPath(std::string p) const {
     struct stat info;
     bool status=true;
     if(p!="")
       {
         if( stat(p.c_str(), &info) != 0 )
         {
             status=false;
             throw std::runtime_error("cannot access "+ p);
         }
         else if( !(info.st_mode & S_IFDIR) )
         {
             status=false;
             throw std::runtime_error("is not a directory: " +p);
         }
     } else {
       std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
       return false;
     }
     return status;
    }
  public:
    std::string golemfit_path =        std::getenv("SNOTPATH") ? std::getenv("SNOTPATH") : "../";
    std::string compact_file_path =        golemfit_path+"/compact_data/";
    std::string squids_files_path =        golemfit_path+"/resources/Fluxes/conventional/";
    std::string prompt_squids_files_path = golemfit_path+"/resources/Fluxes/prompt/";
    std::string neutrino_cc_xs_spline_path= golemfit_path+"/resources/CrossSections/sigma_nu_CC_iso.fits";
    std::string antineutrino_cc_xs_spline_path = golemfit_path+"/resources/CrossSections/sigma_nubar_CC_iso.fits";
    std::string neutrino_nc_xs_spline_path = golemfit_path+"/resources/CrossSections/sigma_nu_NC_iso.fits";
    std::string antineutrino_nc_xs_spline_path = golemfit_path+"/resources/CrossSections/sigma_nubar_NC_iso.fits";
    std::string diff_neutrino_cc_xs_spline_path = golemfit_path+"/resources/CrossSections/dsdxdy_nu_CC_iso.fits";
    std::string diff_antineutrino_cc_xs_spline_path = golemfit_path+"/resources/CrossSections/dsdxdy_nubar_CC_iso.fits";
    std::string diff_neutrino_nc_xs_spline_path = golemfit_path+"/resources/CrossSections/dsdxdy_nu_NC_iso.fits";
    std::string diff_antineutrino_nc_xs_spline_path = golemfit_path+"/resources/CrossSections/dsdxdy_nubar_NC_iso.fits";
    std::string nufate_xs_table_location =  golemfit_path+"/resources/CrossSections/NuFATECrossSections.h5";
    std::string data_path =                golemfit_path+"/data/";
    std::string mc_path =                  golemfit_path+"/monte_carlo/";
    std::string oversize_path =            golemfit_path+"/resources/OversizeCorrections/";
    std::string domeff_spline_path =       golemfit_path+"/resources/DOMEffSplines/";
    std::string hqdomeff_spline_path =       golemfit_path+"/resources/RELDOMEffSplines/";
    std::string initial_flux_files_path =  golemfit_path+"/resources/InitialFluxes/";
    std::string dark_spline_path =  golemfit_path+"/resources/DarkMatterSplines/";
    std::string astrophysical_models_spline_path =  golemfit_path+"/resources/AstrophysicalFluxesSplines/";
    std::string skymaps_spline_path =  golemfit_path+"/resources/FermiSkyMap/";
    std::string holeice_spline_path =  golemfit_path+"/resources/HoleIceSplines/";
    std::string attenuation_spline_path =  golemfit_path+"/resources/AttenuationSplines/";
    std::string anisotropy_spline_path =  golemfit_path+"/resources/TauSplines/";
    std::string conventional_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/conventional_atmospheric.hdf5";
    std::string kaon_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/kaon_atmospheric.hdf5";
    std::string pion_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/pion_atmospheric.hdf5";
    std::string prompt_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/prompt/prompt_atmospheric.hdf5";
    std::string astro_nusquids_file =  golemfit_path+"/resources/Fluxes/astro/astro.hdf5";
    std::string selfveto_spline_path =  golemfit_path+"/resources/SelfVetoSplines/";
    std::string barr_resources_location =  golemfit_path+"/resources/barrgradients/";
    std::string ice_covariance_location =  golemfit_path+"/resources/IceCovarianceMatrix/";
    std::string ice_gradient_location =  golemfit_path+"/resources/IceGradients/";
    std::string atmospheric_density_spline_location =  golemfit_path+"/resources/AtmosphericZenithVariationSplines/";
    std::string atmospheric_kaonlosses_spline_location =  golemfit_path+"/resources/AtmosphericKaonLossesSplines/kaon_loses_1s.fits";
    DataPaths(){};
    DataPaths(std::string golemfit_base_path){SetAllPathsToDefaults(golemfit_base_path);}

    bool CheckDataPaths() const {
      CheckDataPath(compact_file_path);
      CheckDataPath(squids_files_path);
      CheckDataPath(prompt_squids_files_path);
      CheckDataPath(neutrino_cc_xs_spline_path);
      CheckDataPath(antineutrino_cc_xs_spline_path);
      CheckDataPath(neutrino_nc_xs_spline_path);
      CheckDataPath(antineutrino_nc_xs_spline_path);
      CheckDataPath(data_path);
      CheckDataPath(mc_path);
      CheckDataPath(oversize_path);
      CheckDataPath(domeff_spline_path);
      CheckDataPath(hqdomeff_spline_path);
      CheckDataPath(initial_flux_files_path);
      CheckDataPath(nufate_xs_table_location);
      return true;
    }

    void SetAllPathsToDefaults(std::string golemfit_base_path){
       golemfit_path = golemfit_base_path;
       compact_file_path =        golemfit_path+"/compact_data/";
       squids_files_path =        golemfit_path+"/resources/Fluxes/conventional/";
       prompt_squids_files_path = golemfit_path+"/resources/Fluxes/prompt/";
       neutrino_cc_xs_spline_path= golemfit_path+"/resources/CrossSections/sigma_nu_CC_iso.fits";
       antineutrino_cc_xs_spline_path= golemfit_path+"/resources/CrossSections/sigma_nubar_CC_iso.fits";
       neutrino_nc_xs_spline_path= golemfit_path+"/resources/CrossSections/sigma_nu_NC_iso.fits";
       antineutrino_nc_xs_spline_path= golemfit_path+"/resources/CrossSections/sigma_nubar_NC_iso.fits";
       diff_neutrino_cc_xs_spline_path= golemfit_path+"/resources/CrossSections/dsdxdy_nu_CC_iso.fits";
       diff_antineutrino_cc_xs_spline_path= golemfit_path+"/resources/CrossSections/dsdxdy_nubar_CC_iso.fits";
       diff_neutrino_nc_xs_spline_path= golemfit_path+"/resources/CrossSections/dsdxdy_nu_NC_iso.fits";
       diff_antineutrino_nc_xs_spline_path= golemfit_path+"/resources/CrossSections/dsdxdy_nubar_NC_iso.fits";
       nufate_xs_table_location =  golemfit_path+"/resources/CrossSections/NuFATECrossSections.h5";
       data_path =                golemfit_path+"/data/";
       mc_path =                  golemfit_path+"/monte_carlo/";
       oversize_path =            golemfit_path+"/resources/OversizeCorrections/";
       domeff_spline_path =       golemfit_path+"/resources/DOMEffSplines/";
       hqdomeff_spline_path =       golemfit_path+"/resources/RELDOMEffSplines/";
       initial_flux_files_path =  golemfit_path+"/resources/InitialFluxes/";
       dark_spline_path =  golemfit_path+"/resources/DarkMatterSplines/";
       astrophysical_models_spline_path = golemfit_path+"/resources/AstrophysicalFluxesSplines/";
       skymaps_spline_path =  golemfit_path+"/resources/FermiSkyMap/";
       holeice_spline_path =  golemfit_path+"/resources/HoleIceSplines/";
       attenuation_spline_path =  golemfit_path+"/resources/AttenuationSplines/";
       anisotropy_spline_path =  golemfit_path+"/resources/TauSplines/";
       conventional_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/conventional_atmospheric.hdf5";
       kaon_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/kaon_atmospheric.hdf5";
       pion_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/conventional/pion_atmospheric.hdf5";
       prompt_nusquids_atmospheric_file =  golemfit_path+"/resources/Fluxes/prompt/prompt_atmospheric.hdf5";
       astro_nusquids_file =  golemfit_path+"/resources/Fluxes/astro/astro.hdf5";
       selfveto_spline_path =  golemfit_path+"/resources/SelfVetoSplines/";
       barr_resources_location =  golemfit_path+"/resources/barrgradients/";
       ice_covariance_location =  golemfit_path+"/resources/IceCovarianceMatrix/";
       ice_gradient_location =  golemfit_path+"/resources/IceGradients/";
       atmospheric_density_spline_location =  golemfit_path+"/resources/AtmosphericZenithVariationSplines/";
       atmospheric_kaonlosses_spline_location =  golemfit_path+"/resources/AtmosphericKaonLossesSplines/kaon_loses_1s.fits";
    }
};

unsigned int null_categorizer(const Event& e);

unsigned int symmetric_categorizer(const Event& e);

struct SteeringParams {
  bool fastmode=false;
  bool frequentist=true;
  bool add_contribution_from_tau_neutrinos=false;
  double fastmode_scaling = 0.5;
  bool quiet=false;
  bool debug_verbose=false;
  bool paranoid=false;
  double minFitEnergy = pow(10.,4.78); // ~60 TeV
  double maxFitEnergy = 1.e7; // 10 PeV
  double minCosth = -1.;
  double maxCosth = 1;
  double minAzimuth = 0.;
  double maxAzimuth = 2.*boost::math::constants::pi<double>();
  double minLength = 10.;
  double maxLength = 1000.0;
  double logEbinEdge = 4.78;
  double logEbinWidth = 0.111;//0.04530612;
  double cosThbinEdge = -1.0;
  double cosThbinWidth = 0.2;
  double azimuthbinEdge = 0.0;
  double azimuthbinWidth = 2.*boost::math::constants::pi<double>();
  double lengthbinEdge = log10(10.);
  double lengthbinWidth = 0.1;//0.51168558; // previous width was too wide
  DiffuseFitType diffuse_fit_type = DiffuseFitType::SinglePowerLaw;
  double baseline_astro_normalization = 1e-18;
  double baseline_astro_spectral_index = 0;
  double baseline_galactic_astro_normalization = 1e-18;
  double baseline_galactic_astro_spectral_index = 0;
  bool spline_dom_efficiency = true;
  bool spline_hqdom_efficiency = false;
  bool spline_hole_ice = true;
  bool spline_attenuation = false;
  bool spline_anisotropy = false;
  bool apply_anisotropy_systematic = false;
  bool use_nusquids_height_sampling = false;
  bool use_single_conventional_atmospheric_neutrino_file = true;
  double singleCascadeLengthSystematicAmplitude1 = 16.10;
  double singleCascadeLengthSystematicAmplitude2 = -0.82;
  double singleCascadeLengthSystematicOffset1 = 4.91;
  double singleCascadeLengthSystematicOffset2 = -1.68;
  double trackLengthSystematicAmplitude1 = 0;
  double trackLengthSystematicAmplitude2 = 0;
  double trackLengthSystematicOffset1 = 0;
  double trackLengthSystematicOffset2 = 0;
  double doubleCascadeLengthSystematicAmplitude1 = 3.74;
  double doubleCascadeLengthSystematicAmplitude2 = 0.01;
  double doubleCascadeLengthSystematicOffset1 = 1.02;
  double doubleCascadeLengthSystematicOffset2 = -0.15;
  double lengthSystematicAnisotropyPhase = -0.15707963267948966;
  MCType simType {MCType::NuGen};
  sampleTag sampleToLoad {sampleTag::Sterile};
  //sampleTag sampleToLoad {sampleTag::MagicTau};
  size_t evalThreads=4;
  bool use_precalculated_nusquids_fluxes= false;
  std::string atmospheric_modelName = "honda2006";
  std::string atmospheric_kneecorrection_modelName = "gaisserH3a_elbert";
  std::string cosmic_muon_model="GaisserH4a";
  std::string prompt_modelName="BERSS_H3a_central"; //"sarcevic_std";
  bool correct_prompt_knee = true; // in original HESE prompt is not corrected // this should be true
  bool readCompact=false;
  CompactMode compactMode = CompactMode::dump;
  std::vector<unsigned int> years={0,1,2,3,4,5};
  std::map<unsigned int, double> fullLivetime = std::map<unsigned int,double>{{0,331.*86400.}, {1,331.*86400.}, {2,326.385*86400.}, {3,358.4*86400.}, {4,368.308*86400.}, {5,362.525*86400.}, {999,100*365.*86400.}}; //999 for testing long lifetime
  std::vector<double> initial_flavor_ratio {1,1,1};// usual notation
  double fundamental_muongun_scaling=7.409614804931514;//6.348294179520329; // this is to match the muon normalization
  bool do_HESE_reshuffle=true;
  bool allow_mismatched_parameters = false;
  double track_to_shower_missid=0.10; // 10% of the tracks become showers. Magic!
  std::function<unsigned int(const Event&)> categorizer=symmetric_categorizer;
  // for LW simulations only
  std::string simToLoad="p2_0";
  std::string oversizeFunction="NullCorrection";
  std::string pi0_skymap = "fermi_pi0_new.hdf5";
  bool use_legacy_selfveto_calculation = true;
  std::string selfveto_model = "HGH2012_SIBYLL2_3c";
  std::string sterile_model_label = "0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000";
  std::string ice_covariance_filename = "RoughCovarianceMatrix.csv";
  std::vector<std::string> ice_gradient_filename = {"Energy_Eff_Grad_1.txt","Energy_Eff_Grad_2.txt"};
  std::string atmospheric_density_spline_filename = "atm_density_1s.fits";
  bool use_ice_covariance = false;
  bool use_ice_gradients = false;
  bool nusquids_progress_bar = false;
  bool load_barr_gradients = true;
  bool load_atmospheric_density_spline = true;
  bool load_kaon_losses_spline = true;
  bool use_2d_prior = false;
  bool fastmode_detector_systematics_reset = true;
  EarthPropagationInterfaces earth_propagation_interface = EarthPropagationInterfaces::nuSQuIDS;
  AstrophysicalNeutrinoModel ad_hoc_astro_model = AstrophysicalNeutrinoModel::SteckerAGN;
  std::vector<AstrophysicalNeutrinoModel> ad_hoc_astro_models_to_load_ {AstrophysicalNeutrinoModel::SteckerAGN,
                                                                        AstrophysicalNeutrinoModel::Fang_UHECR,
                                                                        AstrophysicalNeutrinoModel::LLGAGN_modelB1,
                                                                        AstrophysicalNeutrinoModel::LLGAGN_modelB4,
                                                                        AstrophysicalNeutrinoModel::LLGAGN_twocompModel,
                                                                        AstrophysicalNeutrinoModel::Maria_BLLacs,
                                                                        AstrophysicalNeutrinoModel::Murase_chockedJets,
                                                                        AstrophysicalNeutrinoModel::SBG_minBmodel,
                                                                        AstrophysicalNeutrinoModel::Tavecchi_lowPower,
                                                                        AstrophysicalNeutrinoModel::TDE_WinterBiehl,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_MixComp_beta2,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_MixComp_beta25 ,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_Proton,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_MixComp_beta25,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_MixComp_beta2,
                                                                        AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_Proton,
                                                                        AstrophysicalNeutrinoModel::Ahlers_GZK};
  std::function<bool(const Event &)> is_sideband = [](const Event & e){return false;};
  SteeringParams(){};
  SteeringParams(sampleTag sample):sampleToLoad(sample){
    if(sample == sampleTag::Sterile){
      simType = MCType::LeptonInjector;
      years = {2011};
      simToLoad = "Platinum_Split";
      //simToLoad = "Platinum_Central_Split";
      logEbinEdge = log10(400);
      logEbinWidth = 0.06;
      minFitEnergy = 4.0e2;
      maxFitEnergy = 2.0e4;
      minCosth = -1;
      maxCosth = 0.0;
      cosThbinEdge = 0.;
      cosThbinWidth = 0.1;
      fastmode = false;
      spline_dom_efficiency = true;
      spline_hqdom_efficiency = true;
      spline_hole_ice = true;
      load_barr_gradients = true;
      load_atmospheric_density_spline = true;
      load_kaon_losses_spline = true;
      use_ice_gradients = true;
      use_precalculated_nusquids_fluxes = true;
      fullLivetime = std::map<unsigned int, double> { {0, 28748982.0}, {2011, 28748982.0}}; // 2011 lifetime
    }
    else if(sample == sampleTag::HESE){
      simToLoad="p2_0";
    }
    else if(sample == sampleTag::MagicTau){
      simToLoad="p2_0";
    }
  }
  void SetSquareSideband(double Emin, double Emax, double costhmin, double costhmax){
    is_sideband = [=](const Event & e){return e.energy > Emin && e.energy < Emax && cos(e.zenith) > costhmin && cos(e.zenith) < costhmax;};
  }
};

struct NewPhysicsParams {
  NewPhysicsType type = NewPhysicsType::None;
  // properties for lorentz violation
  unsigned int index = 0;
  double dm41sq = 0.;
  double th12 = 0;
  double th13 = 0;
  double th23 = 0;
  double th14 = 0;
  double th24 = 0;
  double th34 = 0;
  double del13 = 0;
  double del14 = 0;
  double del24 = 0;
  double lambda_1 = 1.e60;
  double lambda_2 = 1.e60;
  double n_lv = 1;// c-term
  // properties for simple NSI
  double epsilon_mutau = 0;
  double epsilon_prime = 0;
  // properties for dark properties
  double dark_matter_mass = 0;
  double dark_matter_lifetime = std::numeric_limits<double>::max();
  // dark matter interactions and decay channels
  double dark_matter_annihilation_cross_section= 0;
  double dark_matter_neutrino_scattering_mediator_mass = std::numeric_limits<double>::max();
  double dark_matter_neutrino_scattering_coupling = 0;
  std::string dark_matter_annihilation_channel = "anihilation_bbbar_bur_";
  std::string dark_matter_decay_channel = "decay_bbbar_ein_";
  DMnu::Interaction dark_matter_scattering_process = DMnu::Interaction::FV;
  // properties for xs
  double xs0 = 1.;
  double xs1 = 1.;
  double xs2 = 1.;
  double xs3 = 1.;
  /* Probably will move all of this to several class structures and subclasses. TODO CAD. // keep the struct simple and dirty for now.
  // constructor for LV
  NewPhysicsParams(double th12,double th13, double th23, double del13, double lambda):type(NewPhysicsType::LorentzViolation),
    th12(th12),th13(th13),th23(th23),del13(del13),lambda(lambda){}
  // constructor for NSI
  NewPhysicsParams(double epsilon_mutau, double epsilon_prime):type(NewPhysicsType::NonStandardInteraction),
    th12(th12),th13(th13),th23(th23),del13(del13),lambda(lambda){}
  // constructor for dark matter decay
  NewPhysicsParams(double dark_matter_mass, double dark_matter_lifetime, std::string dark_matter_decay_channel):type(NewPhysicsType::DarkMatterDecay),
    dark_matter_mass(dark_matter_mass),dark_matter_lifetime(dark_matter_lifetime),dark_matter_decay_channel(dark_matter_decay_channel){}
  // constructor for dark matter annihilation
  NewPhysicsParams(double dark_matter_mass, double dark_matter_annihilation_cross_section, std::string dark_matter_annihilation_channel):type(NewPhysicsType::DarkMatterAnnihilation),
    dark_matter_mass(dark_matter_mass),dark_matter_annihilation_cross_section(dark_matter_annihilation_cross_section),dark_matter_annihilation_channel(dark_matter_annihilation_channel){}
  // constructor for dark matter scattering
  NewPhysicsParams(double dark_matter_mass, double dark_matter_neutrino_scattering_mediator_mass, double dark_matter_neutrino_scattering_coupling, std::string dark_matter_scattering_process):type(NewPhysicsType::DarkMatterScattering),
    dark_matter_mass(dark_matter_mass),dark_matter_annihilation_cross_section(dark_matter_annihilation_cross_section),dark_matter_annihilation_channel(dark_matter_annihilation_channel){}
  */
  NewPhysicsParams(){};
};

} // namespace golemfit

#endif /* GOLEMPARAMETERS_H_ */
