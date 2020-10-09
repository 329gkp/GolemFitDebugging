#include <stdio.h>
#include <iostream>
#include "GolemFit.h"
#include <nuSQuIDS/marray.h>
#include <fenv.h>

int main(int argc, char* argv[]){
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    using namespace golemfit;

    // Initialize parameter sets at defaults
    std::string base_path = "/home/carguelles/vault/sterilize_cvmfs/sources/GolemFit//";
    DataPaths        dp(base_path);
    //DataPaths        dp("../");
    //dp.barr_resources_location = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes_Final/IceCubeHondaGaisser_th34_0.000_0.0/";
    //dp.barr_resources_location = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Update_H2_HG_SB23C_th34_0.000/";
    //dp.barr_resources_location = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/H2_HG_SB23C_th34_0.000/";
    // THESE two lines are spencer-validates files
    //dp.barr_resources_location = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Flux_AIRS_sib_HG_th24_dm2/";
    //dp.conventional_nusquids_atmospheric_file = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Flux_AIRS_sib_HG_th24_dm2/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5";
    //dp.barr_resources_location = base_path + "/resources/Fluxes/";
    //dp.mc_path = base_path +"/local_monte_carlo";
    //dp.mc_path = "../local_monte_carlo";
    dp.compact_file_path  = "/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/FastMC/Platinum_Split_2.00_7.6_Flux_AIRS_sib_HG_th24_dm2/0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000/";
    SteeringParams   sp(sampleTag::Sterile);
    FitParameters    fp(sampleTag::Sterile);
    //FitParameters    fp(sampleTag::Sterile,FluxComponent::diffuseAstro_mu);
    //FitParameters    fp(sampleTag::Sterile,FluxComponent::atmPrompt);
    FitParameters    fp2(sampleTag::Sterile);
    NewPhysicsParams npp;

    sp.fastmode_scaling = 1;
    sp.fastmode = true;
    //sp.fastmode = false;
    //sp.compactMode = golemfit::CompactMode::jason;
    //sp.readCompact = true;
    sp.use_single_conventional_atmospheric_neutrino_file = true;
    sp.use_ice_gradients = false;
    sp.spline_attenuation = true;
    sp.load_atmospheric_density_spline = true;
    sp.fullLivetime = {{0,(7.*365*24*60*60.)}};
    //sp.is_sideband = [](const Event & e){return (cos(e.zenith) < -0.5);};
    //sp.SetSquareSideband(0.0,0.0,0.0,0.0);

    sp.diffuse_fit_type = DiffuseFitType::AdHocModel;
    sp.ad_hoc_astro_model = AstrophysicalNeutrinoModel::Ahlers_GZK;

    bool fancy_systematics = false;

    sp.spline_dom_efficiency = fancy_systematics;
    sp.spline_hqdom_efficiency = fancy_systematics;
    sp.spline_hole_ice = fancy_systematics;
    sp.load_atmospheric_density_spline = fancy_systematics;
    sp.use_ice_gradients = fancy_systematics;
    sp.load_barr_gradients = fancy_systematics;

    if(sp.fastmode)
      std::cout << "============== Turning on fast mode =============" << std::endl;

    //sp.spline_dom_efficiency = true;
    //sp.spline_hole_ice = true;

    sp.use_nusquids_height_sampling = true;
    //sp.use_nusquids_height_sampling = false;
    sp.use_precalculated_nusquids_fluxes = true;
    //sp.use_precalculated_nusquids_fluxes = false;
    sp.simToLoad = "Platinum_lite";
    //sp.simToLoad = "Platinum_Split";
    //sp.simToLoad = "ares_platinum_welter_split_HE";
    sp.load_barr_gradients = true;

    // Construct object
    std::shared_ptr<GolemFit> golem = std::make_shared<GolemFit>(dp, sp, npp);
    //golem->WriteCompact();
    //while(true) {};
    golem->ReportStatus();
    std::cout << std::endl;

    fp.convNorm = 0.0;
    fp.promptNorm = 0.0;
    fp.astroNorm = 0.0;
    fp.astroNormSec = 1.;
    auto mc_events = golem->SpitMonteCarlo();
    auto weighter = golem->GetEventWeighter(fp);
    double x = 0;
    for(auto event : mc_events){
      double w = weighter(event);
      if(isnan(w))
        std::cout << event << std::endl;
      //std::cout << w << std::endl;
      x += w;
    }
    std::cout << "suma: " << x << std::endl;

    //FitParametersFlag flags(sampleTag::Sterile);
    FitParametersFlag flags(true);
    //fp.astroDeltaGamma = 0.2;
    flags.astroDeltaGamma = false;
    //flags.convNorm = false;
    //flags.barrWP = false;
    //flags.barrWM = false;
    //flags.barrWP=false;
    //flags.zenithCorrection = false;
    golem->SetFitParametersFlag(flags);

    /*
    fp.convNorm = 1.2;
    fp.barrWP = 0.3;
    fp.barrWM = 0.1;
    fp.barrZP = -0.1;
    fp.barrZM = 0.1;
    fp.barrYP = -0.1;
    fp.barrYM = 0.1;
    fp.barrHP = 0.1;
    fp.barrHM = -0.1;
    fp.zenithCorrection = 0.5;
    */

    FitParameters fp_default(sampleTag::Sterile);

    //golem->SpitWeightedMonteCarloToHDF5(fp_default,"meows_heavy_mc_set.hdf5");

    golem->SetupAsimov(fp);

    std::cout << "============== begin make realization =============" << std::endl;
    golem->SpitRealization(fp,0);
    std::cout << "============== end make realization =============" << std::endl;

    nusquids::marray<double,5> expect = golem->GetExpectation(fp);
    double total_number_of_expected_events = std::accumulate(expect.begin(),expect.end(),0.);
    std::cout << "============== RATE ALERT =============" << std::endl;
    std::cout << "Expecting this many events: " << std::setprecision(15) << total_number_of_expected_events << std::endl;
    std::cout << "============== RATE ALERT =============" << std::endl;

    auto expect_all = golem->GetExpectationAll(fp);
    double total_number_of_expected_events_main = std::accumulate(std::get<0>(expect_all).begin(),std::get<0>(expect_all).end(),0.);
    double total_number_of_expected_events_sideband = std::accumulate(std::get<1>(expect_all).begin(),std::get<1>(expect_all).end(),0.);
    std::cout << "============== RATE ALERT SIDEBAND =============" << std::endl;
    std::cout << "Expecting this many events: " << std::setprecision(15) << total_number_of_expected_events_main << " " << total_number_of_expected_events_sideband<< std::endl;
    std::cout << "============== RATE ALERT SIDEBAND =============" << std::endl;

    {
      fp.kaonLosses = 0.3;
      nusquids::marray<double,5> b = golem->GetExpectation(fp);
      double a = std::accumulate(b.begin(),b.end(),0.);
      std::cout << "============== RATE ALERT =============" << std::endl;
      std::cout << "Expecting this many events: " << std::setprecision(15) << a << std::endl;
      std::cout << "============== RATE ALERT =============" << std::endl;
    }

    FitParameters fp_seed1(sampleTag::Sterile);
    fp_seed1.convNorm = 0.9;
    FitParameters fp_seed2(sampleTag::Sterile);
    fp_seed2.convNorm = 1.1;

    golem->SetFitParametersSeed(fp);
    //golem->SetFitParametersSeed({fp_seed1,fp_seed2});

    Priors priors;
    priors.astro1Comp2DCorrelation = 0.0;

    golem->SetFitPriors(priors);

    // Obtain best fit point
    // golem->SetFitParametersSeed(fp);
    auto bfp = golem->MinLLH();

    std::cout << "Obtained this parameters" << std::endl;
    std::cout << bfp.params << std::endl;;
    std::cout << "Obtained this likelihood" << std::endl;
    std::cout << bfp.likelihood << std::endl;;
    std::cout << "Injected point likelihood" << std::endl;
    std::cout << golem->EvalLLH(fp) << std::endl;;

    auto extents = expect.get_extents();
    for(auto it=extents.begin(); it<extents.end(); it++) {
        std::cout << *it << ' ';
    }
    std::cout << std::endl;

    return 0;
}
