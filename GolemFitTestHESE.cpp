#include <stdio.h>
#include <iostream>
#include "GolemFit.h"
#include "compactIO.h"
#include <nuSQuIDS/marray.h>
#include "json.hpp"

int main(int argc, char* argv[]){
    using namespace golemfit;
    // Initialize parameter sets at defaults
    DataPaths        dp;
    SteeringParams   sp;
    NewPhysicsParams npp;
    FitParameters    fp_asimov;

    sp.spline_dom_efficiency = false;
    sp.spline_hole_ice = false;

    sp.sampleToLoad = sampleTag::HESE;
    //sp.sampleToLoad = sampleTag::MagicTau;
    if(sp.sampleToLoad == sampleTag::HESE)
      sp.simToLoad = "p2_0";
    else
      sp.simToLoad = "new";
    //sp.paranoid = true;
    int c = 0;
    while ((c = getopt (argc, argv, "xdf")) != -1)
      switch (c) {
      case 'x':
        // cross section testing
        std::cout << "Doing xs test" << std::endl;
        sp.sampleToLoad = sampleTag::HESE;
        sp.do_HESE_reshuffle = false;
        sp.debug_verbose = true;
        sp.baseline_astro_spectral_index = -2.;
        npp.type = NewPhysicsType::NuSQuIDS;
        npp.xs0 = 1.0;
        npp.xs1 = 1.0;
        npp.xs2 = 1.0;
        npp.xs3 = 5.0;
        break;
      case 'd':
        std::cout << "Doing dark matter test" << std::endl;
        npp.type = NewPhysicsType::DarkMatterDecay;
        npp.dark_matter_mass = 6;

        npp.type = NewPhysicsType::DarkMatterScattering;
        npp.dark_matter_mass = 0.01;
        npp.dark_matter_neutrino_scattering_coupling = 1;
        npp.dark_matter_neutrino_scattering_mediator_mass = 0.01;
        npp.dark_matter_scattering_process = DMnu::Interaction::FV;
        break;
      case 'f':
        std::cout << "Turning on fast mode" << std::endl;
        sp.fastmode_scaling = 0.25;
        sp.fastmode = true;
        break;
      }

    // Construct object
    std::shared_ptr<GolemFit> golem = std::make_shared<GolemFit>(dp, sp, npp);

    //golem->SetupAsimov(fp_asimov);
    //golem->ConstructLikelihoodProblem(golem->GetFitPriors(), golem->GetFitParameters(), golem->GetFitParametersFlag());
    golem->ReportStatus();
    /*
       std::shared_ptr<DataPaths> dp_new(new DataPaths());

       golem->ReConfig<false>(dp_new, std::shared_ptr<SteeringParams>(nullptr), std::shared_ptr<NewPhysicsParams>(nullptr));

       std::shared_ptr<SteeringParams> sp_new(new SteeringParams());

       golem->ReConfig<false>(std::shared_ptr<DataPaths>(nullptr), sp_new, std::shared_ptr<NewPhysicsParams>(nullptr));

       std::shared_ptr<NewPhysicsParams> npp_new(new NewPhysicsParams());

       golem->ReConfig<false>(std::shared_ptr<DataPaths>(nullptr), std::shared_ptr<SteeringParams>(nullptr), npp_new);
       */

    // Obtain best fit point
    auto bfp = golem->MinLLH();
    FitParameters fp = bfp.params;

    //std::cout << "Obtained the following -LLH: " + std::to_string(bfp.likelihood) << std::endl;
    //std::cout << "These are the best fit parameters: " << std::endl;
    //std::cout << bfp.params.astroNorm << std::endl;

    std::cout << bfp.params << std::endl;;
    std::cout << bfp.likelihood << std::endl;;

    nusquids::marray<double,5> expect = golem->GetExpectation(fp);
    //auto expect_it = expect.begin();
    auto extents = expect.get_extents();
    for(auto it=extents.begin(); it<extents.end(); it++) {
        std::cout << *it << ' ';
    }
    std::cout << std::endl;


    return 0;
}
