#ifndef _H_GOLEM_FIT_
#define _H_GOLEM_FIT_
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <functional>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>
#include <memory>

#include <nuSQuIDS/nuSQuIDS.h>
#include <nuFATE/nuFATE.h>

#include <LeptonWeighter/Weighter.h>
#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>
#include <LeptonWeighter/NNFluxInterface.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/histogram.h>
#include <PhysTools/bin_types.h>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include "GolemParameters.h"
#include "Event.h"
#include "analysisWeighting.h"
#include "compactIO.h"
#include "oversizeWeight.h"
#include "self_veto.h"
#include "DarkMatterInteractions.h"
#include "Earth.h"
#include "GolemEnumDefinitions.h"
#include "GolemMCSet.h"
#include "GolemTools.h"
#include "NSI.h"
#include "LV.h"

namespace golemfit{

struct simpleLocalDataWeighter{
   double operator()(const Event& e) const{
     return(e.cachedWeight);}
};

struct simpleLocalDataWeighterConstructor{
    template<typename DataType>
    simpleLocalDataWeighter operator()(const std::vector<DataType> params) const{
        return simpleLocalDataWeighter();
    }
};


using namespace nusquids;
using namespace phys_tools::histograms;
using namespace phys_tools::likelihood;
using HistType = std::tuple< histogram<5,entryStoringBin<std::reference_wrapper<const Event>>> >;
using AuxHistType = std::tuple< histogram<5,entryStoringBin<std::reference_wrapper<const Event>>> >;
using PriorIndices = std::tuple<phys_tools::likelihood::parameters<>,phys_tools::likelihood::parameters<17,18>,phys_tools::likelihood::parameters<20,21>>;
using BasicPrior = FixedSizePriorSet<
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               GaussianPrior,GaussianPrior,GaussianPrior,GaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,
                               LimitedGaussianPrior,LimitedGaussianPrior,LimitedGaussianPrior,
                               LimitedGaussianPrior,GaussianPrior,
                               GaussianPrior,GaussianPrior,
                               GaussianPrior,
                               LimitedGaussianPrior,GaussianPrior>;
using CPrior = ArbitraryPriorType<PriorIndices, BasicPrior, Gaussian2DPrior, Gaussian2DPrior>::type;
using LType=LikelihoodProblem<std::reference_wrapper<const Event>, HistType,simpleLocalDataWeighterConstructor,sterile::WeighterMaker,CPrior,SAYLikelihood,27>;
//using LType=LikelihoodProblem<std::reference_wrapper<const Event>, HistType,simpleLocalDataWeighterConstructor,sterile::WeighterMaker,CPrior,poissonLikelihood,27>;
using hist0_marray=marray<double,5>;
using hist1_marray=marray<double,5>;
using hist_marray=std::tuple<hist0_marray, hist1_marray>;

template<typename ContainerType, typename HistType, typename BinnerType>
  void bin(const ContainerType& data, HistType& hist, const BinnerType& binner){
  for(const Event& event : data)
    binner(hist,event);
}

template<typename ContainerType, typename HistType, typename AuxHistType, typename BinnerType, typename SidebandFunc>
  void bin(const ContainerType& data, HistType& hist, AuxHistType& hist_aux, const BinnerType& binner, const SidebandFunc& is_sideband_func){
  for(const Event& event : data)
    binner(hist,hist_aux,event, is_sideband_func(event));
}

class GolemFit {
  private:
    // All the local configuration variables
    SteeringParams  steeringParams_;
    DataPaths       dataPaths_;

    //hypothesis point we fit to
    NewPhysicsParams new_physics_params_;

    // to store events
    std::deque<Event> metaEvents_;
    std::deque<Event> mainSimulation_;
    std::deque<Event> sample_;

    // for fast mode only
    std::deque<Event> auxSimulation_;

    // histograms
    HistType dataHist_; // analysis data
    AuxHistType auxdataHist_;
    HistType simHist_; // analysis MC is
    AuxHistType auxsimHist_;

    // minimizing objects
    std::vector<FitParameters> fitSeed_;
    FitParameters asimovSeed_;
    FitParametersFlag fixedParams_;
    Priors priors_;

    // weighter object
    sterile::WeighterMaker DFWM;
    std::shared_ptr<LW::Flux> fluxKaon_,fluxPion_,fluxConv_,fluxPrompt_,fluxAstro_,fluxAstroGalactic_;
    std::shared_ptr<LW::CrossSectionFromSpline> xsw_;
    std::vector<std::shared_ptr<LW::Generator>> mcw_;
    std::shared_ptr<LW::Weighter> pionFluxWeighter_;
    std::shared_ptr<LW::Weighter> kaonFluxWeighter_;
    std::shared_ptr<LW::Weighter> convFluxWeighter_;
    std::shared_ptr<LW::Weighter> promptFluxWeighter_;
    std::shared_ptr<LW::Weighter> astroFluxWeighter_;
    OversizeWeighter osw_;

    // Skymaps
    Healpix_Map<double> galactic_template_;
    double skymaps_galactic_template_normalization_;
    std::unique_ptr<double[]> skymap_flux_array_;

    // nuSQuIDS objects
    const squids::Const units;
    const unsigned int numneu = 3;
    const double sclup_ = 1.e18;
    std::shared_ptr<tools::ScaledNeutrinoCrossSections> xsv;

    // aux vectors // bad
    std::vector<double> run_number_;
    std::vector<double> event_number_;
    std::vector<double> subevent_number_;

    /*
    // lepton weighter objects for nuSQuIDS
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_pion;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_kaon;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_prompt;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_astro;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_pion_xs;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_kaon_xs;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_prompt_xs;
    std::unique_ptr<LW::nuSQUIDSAtmFlux<>> lw_flux_astro_xs;
    */

    // lepton weighter objects for nuSQuIDS
    // also related nusquids members for bsm physics
    std::map<FluxComponent,std::shared_ptr<LW::nuSQUIDSAtmFlux<>>> lw_nusquids_fluxes_;
    std::map<FluxComponent,std::shared_ptr<LW::nuSQUIDSAtmFlux<>>> lw_nusquids_fluxes_xs_;
    std::map<FluxComponent,std::shared_ptr<LW::nuSQUIDSAtmFlux<nuSQUIDSNSI>>> lw_nusquids_fluxes_nsi_;
    std::map<FluxComponent,std::shared_ptr<LW::nuSQUIDSAtmFlux<nuSQUIDSLV>>> lw_nusquids_fluxes_lv_;

    // selfveto objects
    self_veto::AnalyticPassingFraction<> selfvetoConventionalCalculator_;
    self_veto::AnalyticPassingFraction<> selfvetoPromptCalculator_;

    // barr parameter components
    std::map<BarrParameter,std::shared_ptr<LW::nuSQUIDSAtmFlux<>>> barr_parameters_nusquids_gradients_;

    // ice covariance matrix
    nusquids::marray<double,2> ice_covariance_;
    std::vector<nusquids::marray<double,2>> ice_gradient_;
    std::pair< std::vector<double>, std::vector<std::vector<double>> > ice_gradient_array_pair_;

    // Status flags
    bool xs_weighter_constructed_ = (false);
    bool astrophysical_model_splines_loaded_ = (false);
    bool flux_weighter_constructed_ = (false);
    bool lepton_weighter_constructed_ = (false);
    bool oversize_weighter_constructed_ = (false);
    bool dom_efficiency_splines_constructed_ = (false);
    bool hqdom_efficiency_splines_constructed_ = (false);
    bool holeice_resources_loaded_ = (false);
    bool attenuation_resources_loaded_ = (false);
    bool selfveto_resources_loaded_ = (false);
    bool data_histogram_constructed_ = (false);
    bool simulation_histogram_constructed_ = (false);
    bool simulation_loaded_ = (false);
    bool simulation_HESE_loaded_ = (false);
    bool simulation_sterile_loaded_ = (false);
    bool simulation_magic_tau_loaded_ = (false);
    bool mc_generation_weighter_constructed_ = (false);
    bool data_loaded_ = (false);
    bool data_HESE_loaded_ = (false);
    bool data_sterile_loaded_ = (false);
    bool likelihood_problem_constructed_ = (false);
    bool simulation_initialized_ = (false);
    bool fastmode_constructed_= (false);
    bool apply_new_physics_ = (false);
    bool apply_dark_matter_interactions_= (false);
    bool apply_cross_section_modfications = (false);
    bool apply_nusquid_modfications = (false);
    bool dark_matter_component_set_ = (true);
    bool did_HESE_shuffle = (false);
    bool skymaps_loaded_ = (false);
    bool nusquids_nsi_fluxes_calculated_= (false);
    bool nusquids_lv_fluxes_calculated_= (false);
    bool using_compact_weighting_ = (false);
    bool asimov_setup_ = (false);
    bool barr_resources_loaded_ = (false);
    bool ice_covariance_loaded_ = (false);
    bool ice_gradient_loaded_ = (false);
    bool atmospheric_density_spline_loaded_ = (false);
    bool kaon_losses_spline_loaded_ = (false);

    // DOM efficiency splines
    std::map<std::pair<FluxComponent,Topology>, std::shared_ptr<splinetable<>>> domefficiencySplines_;
    DOMEfficiencySetter<Event> * domEffSetter_;
    double domEffSpline_minValidValue_ = 0.5;
    double domEffSpline_maxValidValue_ = 1.5;

    // HQ DOM efficiency splines
    std::map<std::pair<FluxComponent,Topology>, std::shared_ptr<splinetable<>>> hqdomefficiencySplines_;
    DOMEfficiencySetter<Event> * hqdomEffSetter_;
    double hqdomEffSpline_minValidValue_ = 0.5;
    double hqdomEffSpline_maxValidValue_ = 1.5;

    // Hole ice splines
    std::map<std::pair<FluxComponent,Topology>, std::shared_ptr<splinetable<>>> holeIceSplines_;
    HoleIceSetter<Event> * holeIceSetter_;
    double holeIceSpline_minValidValue_ = -1.5;
    double holeIceSpline_maxValidValue_ =  0.5;

    // Hole ice splines
    std::map<std::pair<FluxComponent,LW::ParticleType>, std::shared_ptr<splinetable<>>> attenuationSplines_;
    double attenuationSpline_minValidValue_ = 0.5;
    double attenuationSpline_maxValidValue_ = 1.5;

    // Selfveto splines
    std::map<std::pair<FluxComponent,LW::ParticleType>, std::shared_ptr<splinetable<>>> selfVetoSplines_;

    // astrophysical model spline
    std::map<AstrophysicalNeutrinoModel,std::shared_ptr<splinetable<>>> astrophysicalNeutrinoModels_;
    std::map<AstrophysicalNeutrinoModel,std::pair<double,double>> astrophysicalNeutrinoModelsSplinesValidityRanges_;

    // atmopsheric density spline
    std::shared_ptr<splinetable<>> atmosphericDensityUncertaintySpline_ = nullptr;

    // kaon losses spline
    std::shared_ptr<splinetable<>> kaonLossesUncertaintySpline_ = nullptr;

    // likehood problem object
    std::shared_ptr<LType> prob_;
    std::shared_ptr<LType> auxprob_;
  public:
    // Constructor
    GolemFit(DataPaths dataPaths, SteeringParams steeringParams, NewPhysicsParams snp = NewPhysicsParams());

    template<bool constructor>
    void ReConfig(std::shared_ptr<DataPaths> dataPaths, std::shared_ptr<SteeringParams> steeringParams, std::shared_ptr<NewPhysicsParams> new_physics_params) {
      bool reset_steering = constructor || (bool(steeringParams) && (&steeringParams_ != steeringParams.get()));
      bool reset_data = constructor || (bool(dataPaths) && (&dataPaths_ != dataPaths.get()));
      bool reset_npp = constructor || (bool(new_physics_params) && (&new_physics_params_ != new_physics_params.get()));

      // default minimizer seed
      fitSeed_ = {FitParameters()};

      if(!constructor && using_compact_weighting_) {
          simulation_initialized_ = false;
          apply_new_physics_ = false;
          apply_dark_matter_interactions_ = false;
          apply_cross_section_modfications = false;
          apply_nusquid_modfications = false;
      }

      if(!steeringParams->quiet){
        std::cout << "reset_steering: " << reset_steering << std::endl;
        std::cout << "reset_data: " << reset_data << std::endl;
        std::cout << "reset_npp: " << reset_npp << std::endl;
      }

      std::shared_ptr<DataPaths> old_data = std::shared_ptr<DataPaths>(reset_data ? new DataPaths(dataPaths_) : nullptr);
      if(reset_data) {
        dataPaths_ = *dataPaths;
      }

      std::shared_ptr<SteeringParams> old_steering = std::shared_ptr<SteeringParams>(reset_steering ? new SteeringParams(steeringParams_) : nullptr);
      if(reset_steering) {
        steeringParams_ = *steeringParams;
      }

      std::shared_ptr<NewPhysicsParams> old_npp = std::shared_ptr<NewPhysicsParams>(reset_npp ? new NewPhysicsParams(new_physics_params_) : nullptr);
      if(reset_npp) {
        new_physics_params_ = *new_physics_params;
      }

      if(reset_data) {
        if(!steeringParams_.quiet) std::cout<<"GolemFit constructor: checking paths" <<std::endl;
        CheckDataPaths(dataPaths_);
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.spline_dom_efficiency){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->domeff_spline_path != dataPaths_.domeff_spline_path)) || (old_steering->spline_dom_efficiency != steeringParams_.spline_dom_efficiency)) {
            if(!steeringParams_.quiet) std::cout<<"Loading DOM Efficiency Splines" <<std::endl;
            LoadDOMEfficiencySplines();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.spline_hqdom_efficiency){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->hqdomeff_spline_path != dataPaths_.hqdomeff_spline_path)) || (old_steering->spline_hqdom_efficiency != steeringParams_.spline_hqdom_efficiency)) {
            if(!steeringParams_.quiet) std::cout<<"Loading HQ DOM Efficiency Splines" <<std::endl;
            LoadHQDOMEfficiencySplines();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.spline_hole_ice){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->holeice_spline_path != dataPaths_.holeice_spline_path)) || (old_steering->spline_hole_ice != steeringParams_.spline_hole_ice)) {
            if(!steeringParams_.quiet) std::cout<<"Loading HoleIce Splines" <<std::endl;
            LoadHoleIceResources();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.spline_attenuation){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->attenuation_spline_path != dataPaths_.attenuation_spline_path)) || (old_steering->spline_attenuation != steeringParams_.spline_attenuation)) {
            if(!steeringParams_.quiet) std::cout<<"Loading Attenuation Splines" <<std::endl;
            LoadAttenuationResources();
          }
        }
      }
      /*
      if(reset_steering || reset_data) {
        if(steeringParams_.spline_anisotropy){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->anisotropy_spline_path != dataPaths_.anisotropy_spline_path)) || (old_steering->spline_anisotropy != steeringParams_.spline_anisotropy)) {
            if(!steeringParams_.quiet) std::cout<<"Loading Splines" <<std::endl;
            LoadAnisotropyResources();
          }
        }
      }
      */

      if(reset_steering || reset_data) {
        if(steeringParams_.load_barr_gradients){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->barr_resources_location != dataPaths_.barr_resources_location)) || (old_steering->load_barr_gradients != steeringParams_.load_barr_gradients)) {
            if(!steeringParams_.quiet) std::cout<<"Loading Barr gradients" <<std::endl;
            LoadBarrUncertaintiesResources();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.use_ice_covariance){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->ice_covariance_location != dataPaths_.ice_covariance_location)) || (old_steering->use_ice_covariance != steeringParams_.use_ice_covariance)) {
            if(!steeringParams_.quiet) std::cout<<"Loading ice covariance" <<std::endl;
            LoadIceCovariantMatrix();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.use_ice_gradients){
          if(constructor || (reset_steering && (old_steering->years != steeringParams_.years)) || (reset_data && (old_data->ice_covariance_location != dataPaths_.ice_covariance_location)) || (old_steering->use_ice_covariance != steeringParams_.use_ice_covariance)) {
            if(!steeringParams_.quiet) std::cout<<"Loading ice gradients" <<std::endl;
            LoadIceGradient();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.load_atmospheric_density_spline){
          if(constructor || reset_steering || (reset_data && (old_data->atmospheric_density_spline_location != dataPaths_.atmospheric_density_spline_location)) || (old_steering->atmospheric_density_spline_filename!= steeringParams_.atmospheric_density_spline_filename)) {
            if(!steeringParams_.quiet) std::cout<<"Loading atmospheriic density splines" <<std::endl;
            LoadAtmosphericDensityUncertaintiesResources();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.load_kaon_losses_spline){
          if(constructor || reset_steering || (reset_data && (old_data->atmospheric_kaonlosses_spline_location != dataPaths_.atmospheric_kaonlosses_spline_location))) {
            if(!steeringParams_.quiet) std::cout<<"Loading kaon losses splines" <<std::endl;
            LoadKaonEnergyLossesUncertaintiesResources();
          }
        }
      }

      if(reset_steering || reset_data) {
        if(steeringParams_.diffuse_fit_type == DiffuseFitType::AdHocModel){
          LoadAstrophysicalModelSplines();
        }
        if(steeringParams_.diffuse_fit_type == DiffuseFitType::Galactic){
          LoadSkymaps();
        }
      }

      if(reset_steering || reset_data) {
        if(!steeringParams_.quiet) std::cout<<"Constructing DiffuseWeightMaker" <<std::endl;
        DFWM = sterile::WeighterMaker(domefficiencySplines_,hqdomefficiencySplines_,holeIceSplines_,attenuationSplines_,astrophysicalNeutrinoModels_,atmosphericDensityUncertaintySpline_,kaonLossesUncertaintySpline_,ice_gradient_array_pair_,steeringParams_, dataPaths_);
      }

      //bool reset_events = constructor;

      if(reset_data || reset_steering) {
        if(steeringParams_.readCompact){
          if(constructor || (reset_data && (dataPaths_.mc_path != old_data->mc_path ||
              dataPaths_.data_path != old_data->data_path ||
              dataPaths_.compact_file_path != old_data->compact_file_path)) ||
              (reset_steering && (steeringParams_.simToLoad != old_steering->simToLoad))) {
            if(!steeringParams_.quiet) std::cout<<"Loading compact data" <<std::endl;
            LoadCompact();
            //reset_events = true;
          }
        } else {
          if(constructor || (reset_steering && (steeringParams_.sampleToLoad != old_steering->sampleToLoad))) {
            if(!steeringParams_.quiet) std::cout<<"Loading data" <<std::endl;
            LoadData();
          }
          if(constructor || (reset_steering && (steeringParams_.sampleToLoad != old_steering->sampleToLoad || steeringParams_.years != old_steering->years || steeringParams_.fullLivetime != old_steering->fullLivetime))) {
            if(!steeringParams_.quiet) std::cout<<"Loading MC" <<std::endl;
            LoadMC();
            //reset_events = true;
          }
        }
      }

      if(reset_steering || reset_data) {
        if(constructor || steeringParams_.atmospheric_modelName != old_steering->atmospheric_modelName || steeringParams_.atmospheric_kneecorrection_modelName != old_steering->atmospheric_kneecorrection_modelName || steeringParams_.prompt_modelName != old_steering->prompt_modelName || steeringParams_.correct_prompt_knee != old_steering->correct_prompt_knee || steeringParams_.atmospheric_kneecorrection_modelName != old_steering->atmospheric_kneecorrection_modelName || steeringParams_.baseline_astro_normalization != old_steering->baseline_astro_normalization || steeringParams_.baseline_astro_spectral_index != old_steering->baseline_astro_spectral_index) {
          if(!steeringParams_.quiet) std::cout<<"Loading Flux weighter" <<std::endl;
          ConstructFluxWeighter();
          if(!steeringParams_.quiet) std::cout<<"Loading nusquids objects" <<std::endl;
          if(!steeringParams_.use_precalculated_nusquids_fluxes)
            // need interactions turned on for nominal cross section weighting
            ConstructNuSQuIDSObjectsForStandardModel(true);
        }
      }

      if(steeringParams_.simType == MCType::LeptonInjector) {
        if(reset_data) {
          if(!steeringParams_.quiet) std::cout<<"Loading XS" <<std::endl;
          ConstructCrossSectionWeighter();
        }
        if(reset_data || reset_steering) {
          if(!steeringParams_.quiet) std::cout<<"Loading MC weighter" <<std::endl;
          ConstructMonteCarloGenerationWeighter();
        }
        if(reset_data || reset_steering) {
          if(!steeringParams_.quiet) std::cout<<"Loading Lepton weighter" <<std::endl;
          ConstructLeptonWeighter();
        }
        if(constructor || (reset_data && dataPaths_.oversize_path != old_data->oversize_path)) {
          if(!steeringParams_.quiet) std::cout<<"Loading Oversize weighter" <<std::endl;
          ConstructOversizeWeighter();
        }
      }


      if(reset_data || reset_steering || reset_npp) {
        if(reset_data || reset_npp || (reset_steering && (steeringParams_.spline_dom_efficiency != old_steering->spline_dom_efficiency || steeringParams_.simToLoad != old_steering->simToLoad || steeringParams_.initial_flavor_ratio[0] != old_steering->initial_flavor_ratio[0] || steeringParams_.initial_flavor_ratio[1] != old_steering->initial_flavor_ratio[1] || steeringParams_.initial_flavor_ratio[2] != old_steering->initial_flavor_ratio[2] || steeringParams_.cosmic_muon_model != old_steering->cosmic_muon_model || steeringParams_.fundamental_muongun_scaling != old_steering->fundamental_muongun_scaling || using_compact_weighting_))) {
          if(!(steeringParams_.readCompact && constructor)) {
            if(!steeringParams_.quiet) std::cout<<"Weighting MC" <<std::endl;
            WeightMC();
            using_compact_weighting_ = false;
          }

          if(!steeringParams_.quiet) std::cout<<"Making data hist" <<std::endl;
          ConstructDataHistogram();
          if(!steeringParams_.quiet) std::cout<<"Making sim hist" <<std::endl;
          ConstructSimulationHistogram();
          if(steeringParams_.fastmode){
            if(!steeringParams_.quiet) std::cout<<"Caressing the simulation into fast mode" <<std::endl;
            ConstructFastMode();
          }
        }
        if(!steeringParams_.quiet) std::cout<<"Constructing likelihood problem with default settings" <<std::endl;
        ConstructLikelihoodProblem(Priors(steeringParams_.sampleToLoad), {FitParameters(steeringParams_.sampleToLoad)},FitParametersFlag(steeringParams_.sampleToLoad));
      }
    }

    // Check that the directories where files are mean to be exist
    bool CheckDataPaths(DataPaths dp) const;

    // Check a directory exists and throw a relevant error otherwise.
    bool CheckDataPath(std::string p) const;

    // Error trap bad file paths
    std::string CheckedFilePath(std::string) const;

    bool WriteCompact() const;

  protected:
    // to check events
    bool IsEventHealthy(Event& e, bool print_reason = true) const;
    void ForceFitSeedSanity(FitParameters& fitSeed);
    void ForceFitSeedSanity();
    // Functions to load and unload data
    void LoadData();
    void LoadSterileData();
    void LoadMC();
    void LoadSterileMC();
    void LoadCompact();
    void LoadFastCompact();
    void LoadSkymaps();
    void LoadHealpixSkymap();
    void LoadSimpleHDF5Skymap();
    void ClearData();
    void ClearSimulation();
    void ConstructFastMode();
    // loading DOM efficiency splines
    void LoadDOMEfficiencySplines();
    void LoadHQDOMEfficiencySplines();
    void LoadHoleIceResources();
    void LoadAttenuationResources();
    void LoadIceCovariantMatrix();
    void LoadIceGradient();
    void LoadAnisotropyResources();
    void LoadAstrophysicalModelSplines();
    void LoadSelfVetoResources();
    void LoadBarrUncertaintiesResources();
    void LoadAtmosphericDensityUncertaintiesResources();
    void LoadKaonEnergyLossesUncertaintiesResources();
    // functions to construct the weighters
    void ConstructCrossSectionWeighter();
    void ConstructFluxWeighter();
    void ConstructMonteCarloGenerationWeighter();
    void ConstructLeptonWeighter();
    void ConstructOversizeWeighter();
    // Function to initialize the MC weights
    void WeightMC();
    void InitializeSimulationWeights();
    void SetGalacticComponent();
    void InitializeLengthBias();
    // functions to construct the histograms of data and simulation
    void ConstructDataHistogram();
    void ConstructSimulationHistogram();
    // functions to construct the likelihood problem
    void ConstructLikelihoodProblem(Priors priors, std::vector<FitParameters> nuisanceSeed, FitParametersFlag fixedParams);
    // Construct nuSQuIDS objects on the fly
    void ConstructNuSQuIDSObjects();
    template<typename BaseType = nuSQUIDS, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
    void ReconfigureNuSQuIDSObjectSettings(nuSQUIDSAtm<BaseType>& nusquids_atm, bool oscillations = true, double abs_error = 1.e-10, double rel_error = 1.e-10) const{
      nusquids_atm.Set_GSL_step(gsl_odeiv2_step_rkf45);
      nusquids_atm.Set_rel_error(rel_error);
      nusquids_atm.Set_abs_error(abs_error);
      nusquids_atm.Set_IncludeOscillations(oscillations);
      nusquids_atm.Set_h(20*units.km);
      nusquids_atm.Set_ProgressBar(steeringParams_.nusquids_progress_bar);
    }
    void ConstructNuSQuIDSObjectsForStandardModel(bool interactions = false);
    void ConstructNuSQuIDSObjectsForNSI(bool interactions = false);
    void ConstructNuSQuIDSObjectsForLV(bool interactions = false);
    void ConstructNuSQuIDSObjectsForCrossSectionModification(bool xs = true); // this boolean is a legacy thing
    std::map<FluxComponent,nusquids::marray<double,4>> ConstructNuSQuIDSFluxComponentArrays(nusquids::marray<double,1> costh_nodes, nusquids::marray<double,1> energy_nodes) const;
    void CalculateEarthTransferFunction();
    void CalculateEarthTransferFunction_nuFATE();
    void CalculateEarthTransferFunction_nuSQuIDS();
    double EvaluateEarthTransferFunction(Event & e, FluxComponent flux_component) const;
    // Get Jordi delta
    double GetZenithCorrectionScale() const;
    // Convert coordinates
    std::tuple<double,double> EquatorialToGalacticCoordinates(double ra, double dec) const;
  public:
    // Converters between human and vector forms
    std::vector<double> ConvertFitParameters(FitParameters p) const;
    std::vector<bool> ConvertFitParametersFlag(FitParametersFlag pf) const;
    FitParameters ConvertVecToFitParameters(std::vector<double> vecns) const;
  public:
    hist_marray GetRealizationAll(std::vector<double> fit_parameters, int seed) const;
    hist0_marray GetRealization(std::vector<double> fit_parameters, int seed) const;
    hist_marray GetExpectationAll(std::vector<double> fit_parameters) const;
    hist0_marray GetExpectation(std::vector<double> fit_parameters) const;
    marray<double,2> SpitRealization(std::vector<double> fit_parameters, int seed) const;
    marray<double,2> SpitExpectation(std::vector<double> fit_parameters) const;
    template<unsigned int hist_type>
    std::vector<double> PullBinEdges(int dim, const  HistType& h) const;
    void SetupAsimov(std::vector<double> FitParameters);
    bool CheckAsimovSanity() const;
    hist_marray GetRawExpectationAll(std::vector<double> fit_parameters) const;
    hist0_marray GetRawExpectation(std::vector<double> fit_parameters) const;
    hist_marray GetSquareExpectationAll(std::vector<double> fit_parameters) const;
    hist0_marray GetSquareExpectation(std::vector<double> fit_parameters) const;
    hist_marray GetMaxExpectationErrorAll(std::vector<double> fit_parameters) const;
    hist0_marray GetMaxExpectationError(std::vector<double> fit_parameters) const;
    hist_marray GetWeightedExpectationAll(std::function<double(const Event &)> f) const;
    hist0_marray GetWeightedExpectation(std::function<double(const Event &)> f) const;
    nusquids::marray<std::vector<Event>,5> GetEventInSimulationHistogram() const;
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyZenithExpectationHistogram(std::vector<double> fit_params, int topology = -1) const;
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyLogLengthExpectationHistogram(std::vector<double> fit_params, int topology = -1) const;
    histogram<1,entryStoringBin<std::reference_wrapper<const Event>>> GetLogLengthExpectationHistogram(std::vector<double> fit_params, int topology = -1) const;
  public:
    // Methods to spit out and swallow event samples
    marray<double,2> SpitData() const;
    std::deque<Event> SpitMonteCarlo() const;
    marray<double,2> SpitRealization(FitParameters nuisance, int seed) const;
    marray<double,2> SpitExpectation(FitParameters nuisance) const;
    void SpitWeightedMonteCarloToHDF5(FitParameters nuisance, std::string filename) const;
    double Swallow(marray<double,2> Data);
    bool SetupAsimov(FitParameters nuisance);
    double SetupDataChallenge(int seed, FitParameters fit_parameters);
    double SetupDataChallenge(int seed, FitParameters fit_parameters, NewPhysicsParams new_physics_params);

    // Methods to get histogram binning
    std::vector<double> GetEnergyBinsData() const;
    std::vector<double> GetEnergyBinsMC() const;
    std::vector<double> GetZenithBinsData() const;
    std::vector<double> GetZenithBinsMC() const;
    std::vector<double> GetAzimuthBinsData() const;
    std::vector<double> GetAzimuthBinsMC() const;
    std::vector<double> GetTopologyBinsData() const;
    std::vector<double> GetTopologyBinsMC() const;
    std::vector<double> GetYearBinsData() const;
    std::vector<double> GetYearBinsMC() const;
    std::vector<std::vector<std::vector<double>>> GetBinsData() const;
    std::vector<std::vector<std::vector<double>>> GetBinsMC() const;
    size_t GetTotalNumberOfBins() const;

    // Methods to change the steering parameters
    void SetEnergyBinning();
    void SetZenithBinning();
    void SetAzimuthBinning();

    // functions to check the status of the object
    bool CheckDataLoaded() const                       {return data_loaded_;};
    bool CheckSimulationLoaded() const                 {return simulation_loaded_;};
    bool CheckDOMEfficiencySplinesConstructed() const  {return dom_efficiency_splines_constructed_;};
    bool CheckCrossSectionWeighterConstructed() const  {return xs_weighter_constructed_;};
    bool CheckFluxWeighterConstructed() const          {return flux_weighter_constructed_;};
    bool CheckOversizeWeighterConstructed() const      {return oversize_weighter_constructed_;};
    bool CheckLeptonWeighterConstructed() const        {return lepton_weighter_constructed_;};
    bool CheckDataHistogramConstructed() const         {return data_histogram_constructed_;};
    bool CheckSimulationHistogramConstructed() const   {return simulation_histogram_constructed_;};
    bool CheckLikelihoodProblemConstruction() const    {return likelihood_problem_constructed_;};
    void ReportStatus() const;

    // functions to obtain distributions
    std::function<double(const Event&)> GetEventWeighter(FitParameters fp) const;
    hist_marray GetDataDistributionAll() const;
    hist0_marray GetDataDistribution() const;
    hist_marray GetExpectationAll(FitParameters fp) const;
    hist0_marray GetExpectation(FitParameters fp) const;
    hist_marray GetSquareExpectationAll(FitParameters fit_params) const;
    hist0_marray GetSquareExpectation(FitParameters fit_params) const;
    hist_marray GetMaxExpectationErrorAll(FitParameters fit_params) const;
    hist0_marray GetMaxExpectationError(FitParameters fit_params) const;
    hist_marray GetRealizationAll(FitParameters fp, int seed) const;
    hist0_marray GetRealization(FitParameters fp, int seed) const;
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyZenithExpectationHistogram(FitParameters fit_params, int topology = -1) const;
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyLogLengthExpectationHistogram(FitParameters fit_params, int topology = -1) const;
    histogram<1,entryStoringBin<std::reference_wrapper<const Event>>> GetLogLengthExpectationHistogram(FitParameters fit_params, int topology = -1) const;
    // functions to evaluate the likelihood
    double EvalLLH(std::vector<double> fp, bool include_prior) const;
    double EvalLLH(std::vector<double> fp) const {return EvalLLH(fp,true);}
    double operator()(std::vector<double> fp) const {return EvalLLH(fp);}
    double EvalLLH(FitParameters fp, bool include_prior) const;
    double EvalLLH(FitParameters fp) const {return EvalLLH(fp,true);}
    double operator()(FitParameters fp) const {return EvalLLH(fp);}
    phys_tools::autodiff::FD<27> EvalLLHGradient(std::vector<phys_tools::autodiff::FD<27>> v) const;
    FitResult MinLLH() const {return MinLLH(false);}
    FitResult MinLLH(bool evaluate_sideband) const;
  private:
    // Nasty template part of fit function
    template<typename LikelihoodType>
      bool DoFitLBFGSB(LikelihoodType& likelihood, phys_tools::lbfgsb::LBFGSB_Driver& minimizer) const{
      using namespace phys_tools::likelihood;
      bool succeeded=minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
      return succeeded;
    }

  public:
    // get/set functions

    SteeringParams       GetSteeringParams()    { return steeringParams_;};
    DataPaths            GetDataPaths()         { return dataPaths_;};
    NewPhysicsParams     GetNewPhysicsParams()  { return new_physics_params_;};
    std::vector<FitParameters>        GetFitParametersSeed() { return fitSeed_;};
    Priors               GetFitPriors()         { return priors_;};
    FitParametersFlag    GetFitParametersFlag() { return fixedParams_;};

    void       SetSteeringParams(SteeringParams p)   { steeringParams_=p; DFWM = sterile::WeighterMaker(domefficiencySplines_,hqdomefficiencySplines_,holeIceSplines_,attenuationSplines_,astrophysicalNeutrinoModels_,atmosphericDensityUncertaintySpline_,kaonLossesUncertaintySpline_,ice_gradient_array_pair_,p, dataPaths_);}
    void       SetDataPaths(DataPaths p)             { dataPaths_=p; CheckDataPaths(p); DFWM = sterile::WeighterMaker(domefficiencySplines_,hqdomefficiencySplines_,holeIceSplines_,attenuationSplines_,astrophysicalNeutrinoModels_,atmosphericDensityUncertaintySpline_,kaonLossesUncertaintySpline_,ice_gradient_array_pair_,steeringParams_,p);}
    void       SetNewPhysicsParams(NewPhysicsParams npp);
    void       SetFitParametersSeed(FitParameters fp) { fitSeed_= std::vector<FitParameters>{fp}; ForceFitSeedSanity(); ConstructLikelihoodProblem(priors_, fitSeed_, fixedParams_);}
    void       SetFitParametersSeed(std::vector<FitParameters> fp) { fitSeed_= fp; ForceFitSeedSanity(); ConstructLikelihoodProblem(priors_, fitSeed_, fixedParams_);}
    void       SetFitParametersFlag(FitParametersFlag fpg) { fixedParams_ = fpg; ConstructLikelihoodProblem(priors_, fitSeed_, fixedParams_);}
    void       SetFitParametersPriors(Priors priors) { priors_ = priors; ConstructLikelihoodProblem(priors_, fitSeed_, fixedParams_);}
    void       SetFitPriors(Priors p) { SetFitParametersPriors(p); }

    double GetLivetime() const {
      double total_livetime = 0;
      for(auto year : steeringParams_.years){
        total_livetime+=steeringParams_.fullLivetime.find(year)->second;
      }
      return total_livetime;
    }

    unsigned int GetNumberOfLoadedMCEvents() const { return mainSimulation_.size();}
    std::string GetNeutrinoMCTag() const { return steeringParams_.simToLoad;}

};

} // close namespace SterileSearch

#endif
