#include <cmath>
#include <algorithm>

#include "GolemFit.h"
#include "FastMode.h"
#include "GolemMCSpecifications.h"

namespace golemfit {

/*************************************************************************************************************
 * Constructor
 * **********************************************************************************************************/

GolemFit::GolemFit(DataPaths dataPaths, SteeringParams steeringParams, NewPhysicsParams new_physics_params){
  ReConfig<true>(std::shared_ptr<DataPaths>(new DataPaths(dataPaths)),
      std::shared_ptr<SteeringParams>(new SteeringParams(steeringParams)),
      std::shared_ptr<NewPhysicsParams>(new NewPhysicsParams(new_physics_params)));
}

/*************************************************************************************************************
 * Implementation auxiliary functions
 * **********************************************************************************************************/

void binner (HistType& h, AuxHistType& h_aux, const Event& e, const bool& is_sideband) {
    if(e.topology == 1 and not is_sideband) {
        assert(e.topology >= 0);
        assert(!std::isnan(e.energy));
        assert(!std::isnan(e.zenith));
        assert(!std::isnan(e.azimuth));
        assert(!std::isnan(e.topology));
        assert(!std::isnan(e.year));
        assert(e.year == 0);
        std::get<0>(h).add(e.energy,cos(e.zenith),e.azimuth,e.topology,e.year,amount(std::cref(e)));
    } else if (e.topology == 1 and is_sideband) {
        assert(e.topology >= 0);
        assert(!std::isnan(e.energy));
        assert(!std::isnan(e.zenith));
        assert(!std::isnan(e.azimuth));
        assert(!std::isnan(e.topology));
        assert(!std::isnan(e.year));
        assert(e.year == 0);
        std::get<0>(h_aux).add(e.energy,cos(e.zenith),e.azimuth,e.topology,e.year,amount(std::cref(e)));
    } else {
      throw std::runtime_error("This should never happen in MEOWS");
    }
};

/*************************************************************************************************************
 * Auxiliary function to debug events
 * **********************************************************************************************************/

bool GolemFit::IsEventHealthy(Event& e, bool print_reason) const {
  if(e.primaryEnergy < 0 or e.primaryEnergy == 0 or std::isnan(e.primaryEnergy)){
    if(print_reason) std::cout << "Unhealthy primary energy" << std::endl;
      return false;
  }
  if(e.energy< 0 or e.energy== 0 or std::isnan(e.energy)){
    if(print_reason) std::cout << "Unhealthy energy" << std::endl;
      return false;
  }
  if(std::isnan(e.primaryZenith)){
    if(print_reason) std::cout << "Unhealthy primary zenith" << std::endl;
      return false;
  }
  if(std::isnan(e.zenith)){
    if(print_reason) std::cout << "Unhealthy zenith" << std::endl;
      return false;
  }
  if(std::isnan(e.primaryAzimuth)){
    if(print_reason) std::cout << "Unhealthy primary azimuth" << std::endl;
      return false;
  }
  if(std::isnan(e.azimuth)){
    if(print_reason) std::cout << "Unhealthy azimuth" << std::endl;
      return false;
  }
  if(std::isnan(e.cachedLivetime) or e.cachedLivetime <0){
    if(print_reason) std::cout << "Unhealthy cachedlifetime" << std::endl;
      return false;
  }
  if(std::isnan(e.cachedConvKaonWeight) or e.cachedConvKaonWeight<0){
    if(print_reason) std::cout << "Unhealthy cachedconvkaonweight" << std::endl;
      return false;
  }
  if(std::isnan(e.cachedConvPionWeight) or e.cachedConvPionWeight<0){
    if(print_reason) std::cout << "Unhealthy cachedconvpionweight" << std::endl;
      return false;
  }
  if(std::isnan(e.cachedPromptWeight) or e.cachedPromptWeight<0){
    if(print_reason) std::cout << "Unhealthy cachedpromptweight" << std::endl;
      return false;
  }
  if(std::isnan(e.cachedWeight) or e.cachedWeight<0){
    if(print_reason) std::cout << "Unhealthy cachedweight" << std::endl;
      return false;
  }
  // properties that are not needed for muon gun simulation
  /*
  if(not (e.primaryType == LW::ParticleType::MuPlus or e.primaryType ==  LW::ParticleType::MuMinus or e.primaryType ==  LW::ParticleType::unknown)){
    if(std::isnan(e.x)){
      if(print_reason) std::cout << "Unhealthy x" << std::endl;
        return false;
    }
    if(std::isnan(e.y)){
      if(print_reason) std::cout << "Unhealthy y" << std::endl;
        return false;
    }
    if(std::isnan(e.z)){
      if(print_reason) std::cout << "Unhealthy z" << std::endl;
        return false;
    }
  }
  */

  return true;
}

/*************************************************************************************************************
 * Functions to read and write data
 * **********************************************************************************************************/

void GolemFit::LoadData(){
  auto sampleToLoad = steeringParams_.sampleToLoad;
  if(sampleToLoad == sampleTag::Sterile){
    LoadSterileData();
  }

  data_loaded_ = data_sterile_loaded_;
  if(not data_loaded_)
    throw std::runtime_error("Sample not found.");
  if(data_loaded_)
      asimov_setup_ = false;
}

void GolemFit::LoadSterileData(){
  using namespace phys_tools::tableio;
  using namespace H5Load::sterile;
  try{
    auto dataAction = [&](RecordID id, Event& e, int dataYear){
      if(e.check(false,Level::neutrino)){
        e.sample=sampleTag::Sterile;
        e.year=0;
        e.cachedWeight=1.;
        e.topology = (unsigned int)Topology::track;// its a track
        sample_.push_back(e);
      }
    };
    auto ic86Action=[&](RecordID id, Event& e){ dataAction(id,e,0); };
    readFile(CheckedFilePath(dataPaths_.data_path+"IC86.h5"),ic86Action);
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
  }
  if(!steeringParams_.quiet)
    std::cout << "Loaded " << sample_.size() << " experimental events" << std::endl;
  data_sterile_loaded_=true;
}

void GolemFit::LoadMC(){
  auto sampleToLoad = steeringParams_.sampleToLoad;
  if(sampleToLoad == sampleTag::Sterile){
    if(!steeringParams_.quiet)
      std::cout << "Begin Loading Sterile MC" << std::endl;
    LoadSterileMC();
    if(!steeringParams_.quiet)
      std::cout << "End Loading Sterile MC" << std::endl;
  }
  simulation_loaded_ = simulation_sterile_loaded_;
  if(not simulation_loaded_)
    throw std::runtime_error("Tag simulation not found.");
}

/*************************************************************************************************************
 * Functions to read and write MC
 * **********************************************************************************************************/

void GolemFit::LoadSterileMC(){
    using namespace phys_tools::tableio;

    if(not dom_efficiency_splines_constructed_ and steeringParams_.spline_dom_efficiency)
      throw std::runtime_error("MC cannot be loaded until dom splines are loaded");
    if(not holeice_resources_loaded_ and steeringParams_.spline_hole_ice)
      throw std::runtime_error("MC cannot be loaded until holeice splines are loaded");

    mainSimulation_.clear();
    metaEvents_.clear();

    run_number_.clear();
    event_number_.clear();
    subevent_number_.clear();

    std::map<unsigned int,double> livetime;
    livetime=steeringParams_.fullLivetime;

    std::vector<MCSet> simSetsToLoad = sterile::GetSimulationSets(steeringParams_.simToLoad, dataPaths_.mc_path+ "/" + GetSampleName(steeringParams_.sampleToLoad));

    try{
      auto simAction=[&](RecordID id, Event& e, const MCSet& mc_set, const DOMEfficiencySetter<Event>& domEff, const HoleIceSetter<Event>& holeIce){
        if(e.check(false,Level::neutrino) && e.energy>1){
          unsigned int simYear=std::dynamic_pointer_cast<LW::RangeGenerator>(mc_set.generators.front())->GetSimulationDetails().Get_Year();
          e.sample=sampleTag::Sterile;
          e.year=0;
          //e.isNuGen = false;
          e.cachedLivetime=livetime.find(simYear)->second;
          e.cachedConvPionWeight=0;
          e.cachedConvKaonWeight=0;
          e.cachedConvWeight=0;
          e.cachedPromptWeight=0;
          e.cachedAstroMuWeight=0;
          e.oneWeight=0;

          e.topology = (unsigned int)Topology::track;// its a track

          run_number_.push_back(id.run);
          event_number_.push_back(id.event);
          subevent_number_.push_back(id.subEvent);

          if(steeringParams_.spline_dom_efficiency){
            domEff.setCache(e);
          }

          if(steeringParams_.spline_hole_ice){
            holeIce.setCache(e);
          }

          e.number_of_generated_mc_events = mc_set.number_of_generated_events;

          if(e.primaryType==LW::ParticleType::NuTau || e.primaryType==LW::ParticleType::NuTauBar){
            assert(e.cachedConvPionWeight==0.0);
            assert(e.cachedConvKaonWeight==0.0);
            assert(e.cachedConvWeight==0.0);
            assert(e.cachedPromptWeight==0.0);
            assert(e.cachedAstroMuWeight==0.0);
          }


          mainSimulation_.push_back(e);
        }
      };

      for(auto simSet : simSetsToLoad){
        if(steeringParams_.spline_dom_efficiency){
          if(domefficiencySplines_.empty())
            throw std::runtime_error("LoadSterileMC: DOM efficiency splines are empty. Cannot set cached efficiency.");
          domEffSetter_=new DOMEfficiencySetter<Event>(domefficiencySplines_,simSet.unshadowedFraction);
        }
        if(steeringParams_.spline_hole_ice){
          if(holeIceSplines_.empty())
            throw std::runtime_error("LoadSterileMC: Hole ice splines are empty. Cannot set cached efficiency.");
          holeIceSetter_=new HoleIceSetter<Event>(holeIceSplines_,simSet.holeiceForward);
        }

        auto callback=[&](RecordID id, Event& e){ simAction(id,e,simSet,*domEffSetter_,*holeIceSetter_); };
        if(not steeringParams_.quiet){
          std::cout << "Number of MC sets and splits: " << std::get<0>(simSet.split) << " " << std::get<1>(simSet.split) << std::endl;
        }

        if(std::get<0>(simSet.split)){
          if(!steeringParams_.quiet)
            std::cout << "Reading split Monte Carlo set" << std::endl;
          for(unsigned int split = 0; split < std::get<1>(simSet.split); split++){
            auto path=CheckedFilePath(simSet.path+"/"+simSet.filename+"/"+simSet.filename+"_"+std::to_string(split)+".h5");
            if(not steeringParams_.quiet){
              std::cout << path << std::endl;
            }
            H5Load::sterile::readFile(path,callback);
          }
        } else {
          if(!steeringParams_.quiet)
            std::cout << "Reading unsplit Monte Carlo set" << std::endl;
          auto path=CheckedFilePath(simSet.path+"/"+simSet.filename);
          H5Load::sterile::readFile(path,callback);
        }
      }
    } catch(std::exception& ex) {
      std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
    }
    if(!steeringParams_.quiet)
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set" << std::endl;
    simulation_sterile_loaded_=true;
    assert(not mainSimulation_.empty());
}

/*************************************************************************************************************
 * Functions to do load compact
 * **********************************************************************************************************/

bool GolemFit::WriteCompact() const {
  // ALEJO: be aware that the sideband MC and data are not dump on the compact format.
  if(steeringParams_.compactMode == CompactMode::dump){
    std::string file_path_compact = dataPaths_.compact_file_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + ".meows";
    golemfit::dump::splatData(file_path_compact,0,sample_,(fastmode_constructed_)?metaEvents_:mainSimulation_);
  } else if (steeringParams_.compactMode == CompactMode::jason){
    std::string file_path_compact = dataPaths_.compact_file_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + ".json";
    golemfit::jason::splatData(file_path_compact,sample_,(fastmode_constructed_)?metaEvents_:mainSimulation_);
  } else {
    throw std::runtime_error("Compact mode requested is not implemented. Logic error.");
  }
  return true;
}

void GolemFit::LoadCompact() {
  // ALEJO: be aware that the sideband MC and data are not dump on the compact format.
  if(steeringParams_.compactMode == CompactMode::dump){
    std::string file_path_compact = dataPaths_.compact_file_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + ".meows";
    golemfit::dump::unsplatData(CheckedFilePath(file_path_compact),0,sample_,mainSimulation_);
  } else if (steeringParams_.compactMode == CompactMode::jason){
    std::string file_path_compact = dataPaths_.compact_file_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + ".json";
    golemfit::jason::unsplatData(file_path_compact,sample_,mainSimulation_);
  } else {
    throw std::runtime_error("Compact mode requested is not implemented. Logic error.");
  }
  data_loaded_ = true;
  simulation_loaded_ = true;
  ConstructDataHistogram();
  ConstructSimulationHistogram();
  using_compact_weighting_ = true;
}

/*
bool GolemFit::WriteBabyCompact() const {
  using namespace golemfit::dump;
  std::deque<babyEvent> main_simulation_baby_eventos(mainSimulation_.size());
  for(auto const event& : mainSimulation_){
    babyEvent.push_back(event.spitbabyEvent());
  }
  std::deque<babyEvent> sample_baby_eventos(sample_.size());
  for(auto const event& : sample_){
    sample_baby_eventos.push_back(event.spitbabyEvent());
  }
  std::string file_path_compact = dataPaths_.compact_file_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + ".baby";
  splatData(file_path_compact,sample_baby_eventos,main_simulation_baby_eventos);
}
*/

void GolemFit::ClearData(){
  sample_.clear();
}

void GolemFit::ClearSimulation(){
  mainSimulation_.clear();
  metaEvents_.clear();
}

/*************************************************************************************************************
 * Functions to load to load astrophysical models splines
 * **********************************************************************************************************/

void GolemFit::LoadAstrophysicalModelSplines(){
  for(auto model : steeringParams_.ad_hoc_astro_models_to_load_){
    astrophysicalNeutrinoModels_.insert({model,
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.astrophysical_models_spline_path+"/"+GetAstrophysicalNeutrinoModelName(model)+".fits")))});
    astrophysicalNeutrinoModelsSplinesValidityRanges_.insert({model,
        {astrophysicalNeutrinoModels_[model]->lower_extent(0),astrophysicalNeutrinoModels_[model]->upper_extent(0)}
        });
    if(not steeringParams_.quiet)
      std::cout << "Loaded ad hoc astrophysical model spline: " << GetAstrophysicalNeutrinoModelName(model) << std::endl;
  }
  assert(astrophysicalNeutrinoModels_.size() != 0);
  astrophysical_model_splines_loaded_ = true;
}

/*
void GolemFit::LoadAstrophysicalModelSplines(){
  for(auto model : { AstrophysicalNeutrinoModel::SteckerAGN }){
    astrophysicalNeutrinoModels_.insert({model,
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.astrophysical_models_spline_path+"/"+GetAstrophysicalNeutrinoModelName(model)+".fits")))});
    astrophysicalNeutrinoModelsSplinesValidityRanges_.insert({model,
        {astrophysicalNeutrinoModels_[model]->lower_extent(0),astrophysicalNeutrinoModels_[model]->upper_extent(0)}
        });
  }
  astrophysical_model_splines_loaded_ = true;
}
*/

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

void GolemFit::LoadDOMEfficiencySplines(){
  domefficiencySplines_.clear();
  if(!steeringParams_.quiet) std::cout << "Loading DOM efficiency splines..." << std::endl;

  auto register_domefficiency_spline = [&](FluxComponent component, Topology topology) {
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
    if(!steeringParams_.quiet)
      std::cout << "Loading " << dataPaths_.domeff_spline_path+"/"+sample_name+"/"+sample_name+"_domefficiency_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits" << std::endl;
    domefficiencySplines_.insert({std::make_pair(component,topology),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.domeff_spline_path+"/"+sample_name+"/"+sample_name+"_domefficiency_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits")))});
  };

  register_domefficiency_spline(FluxComponent::atmConv,Topology::track);
  register_domefficiency_spline(FluxComponent::atmPrompt,Topology::track);
  register_domefficiency_spline(FluxComponent::diffuseAstro_mu,Topology::track);

  if(steeringParams_.sampleToLoad == sampleTag::HESE){
    register_domefficiency_spline(FluxComponent::atmConv,Topology::shower);
    register_domefficiency_spline(FluxComponent::atmPrompt,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_e,Topology::track);
    register_domefficiency_spline(FluxComponent::diffuseAstro_e,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_mu,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_tau,Topology::track);
    register_domefficiency_spline(FluxComponent::diffuseAstro_tau,Topology::shower);
  }

  assert(domefficiencySplines_.size() != 0);

  // assuming all splines have the same validity ranges
  domEffSpline_minValidValue_ = domefficiencySplines_.begin()->second->lower_extent(2);
  domEffSpline_maxValidValue_ = domefficiencySplines_.begin()->second->upper_extent(2);

  if(!steeringParams_.quiet) std::cout << "DOM efficiency splines validity range set to [" + std::to_string(domEffSpline_minValidValue_) + ", " + std::to_string(domEffSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetDOMEfficiencySplines(domefficiencySplines_);

  dom_efficiency_splines_constructed_=true;
}

/*************************************************************************************************************
 * Functions to load to load HQDOM efficiency splines
 * **********************************************************************************************************/

void GolemFit::LoadHQDOMEfficiencySplines(){
  hqdomefficiencySplines_.clear();
  if(!steeringParams_.quiet) std::cout << "Loading HQ DOM efficiency splines..." << std::endl;

  auto register_domefficiency_spline = [&](FluxComponent component, Topology topology) {
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
      hqdomefficiencySplines_.insert({std::make_pair(component,topology),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.hqdomeff_spline_path+"/"+sample_name+"/"+sample_name+"_hqdomefficiency_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits")))});
  };

  register_domefficiency_spline(FluxComponent::atmConv,Topology::track);
  register_domefficiency_spline(FluxComponent::atmPrompt,Topology::track);
  register_domefficiency_spline(FluxComponent::diffuseAstro_mu,Topology::track);

  if(steeringParams_.sampleToLoad == sampleTag::HESE){
    register_domefficiency_spline(FluxComponent::atmConv,Topology::shower);
    register_domefficiency_spline(FluxComponent::atmPrompt,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_e,Topology::track);
    register_domefficiency_spline(FluxComponent::diffuseAstro_e,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_mu,Topology::shower);
    register_domefficiency_spline(FluxComponent::diffuseAstro_tau,Topology::track);
    register_domefficiency_spline(FluxComponent::diffuseAstro_tau,Topology::shower);
  }

  assert(hqdomefficiencySplines_.size() != 0);

  // assuming all splines have the same validity ranges
  hqdomEffSpline_minValidValue_ = hqdomefficiencySplines_.begin()->second->lower_extent(2);
  hqdomEffSpline_maxValidValue_ = hqdomefficiencySplines_.begin()->second->upper_extent(2);

  if(!steeringParams_.quiet) std::cout << "HQ DOM efficiency splines validity range set to [" + std::to_string(hqdomEffSpline_minValidValue_) + ", " + std::to_string(hqdomEffSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetHQDOMEfficiencySplines(hqdomefficiencySplines_);

  hqdom_efficiency_splines_constructed_=true;
}
/*************************************************************************************************************
 * Functions to load to load selfveto splines
 * **********************************************************************************************************/

void GolemFit::LoadSelfVetoResources(){
  selfVetoSplines_.clear();
  if(!steeringParams_.quiet) std::cout << "Loading Self-veto splines..." << std::endl;

  auto register_selfveto_spline = [&](FluxComponent component, LW::ParticleType particle){
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
    selfVetoSplines_.insert({std::make_pair(component,particle),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.selfveto_spline_path+"/"+ steeringParams_.selfveto_model + "_nuveto_"+GetFluxComponentName(component)+"_"+GetParticleName(particle)+".fits")))});
  };

  if(steeringParams_.sampleToLoad == sampleTag::HESE or steeringParams_.sampleToLoad == sampleTag::MagicTau){
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuE);
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuEBar);
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuMu);
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuMuBar);
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuTau); // secretly its nue. love. ca.
    register_selfveto_spline(FluxComponent::atmConv,LW::ParticleType::NuTauBar); // secretly its nuebar. love. ca.
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuE);
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuEBar);
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuMu);
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuMuBar);
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuTau); // secretly its nue. love. ca.
    register_selfveto_spline(FluxComponent::atmPrompt,LW::ParticleType::NuTauBar); // secretly its nuebar. love. ca.
  }

  assert(selfVetoSplines_.size() != 0);
  selfveto_resources_loaded_ = true;
}

void GolemFit::LoadIceCovariantMatrix(){
  if(!steeringParams_.quiet) std::cout << "Loading covariance matrix..." << std::endl;
  //ice_covariance.clear(); // missing feature in marray; need to add CA.
  std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
  ice_covariance_ = nusquids::quickread(CheckedFilePath(dataPaths_.ice_covariance_location + "/" + sample_name + "/" + steeringParams_.ice_covariance_filename));
  assert(not ice_covariance_.empty());
  ice_covariance_loaded_ = true;
}

void GolemFit::LoadIceGradient(){
  if(!steeringParams_.quiet) std::cout << "Loading ice gradients..." << std::endl;
  ice_gradient_.clear();
  std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
  for(std::string ice_gradient_filename : steeringParams_.ice_gradient_filename){
    ice_gradient_.push_back(nusquids::quickread(CheckedFilePath(dataPaths_.ice_gradient_location + "/" + sample_name + "/" + ice_gradient_filename)));
  }
  assert(not ice_gradient_.empty());

  auto fgrad = ice_gradient_.front();
  std::vector<double> energy_bins; // edges
  for(unsigned int i = 0; i < fgrad.extent(0); i++){
    energy_bins.push_back(fgrad[i][0]);
  }
  energy_bins.push_back(fgrad[fgrad.extent(0)-1][1]);
  assert(energy_bins.size() > 0);

  std::vector<std::vector<double>> simple_grads;
  for(unsigned int i = 0; i < ice_gradient_.size(); i++){
    std::vector<double> grad;
    for(unsigned int j = 0; j < ice_gradient_[i].extent(0); j++){
      grad.push_back(ice_gradient_[i][j][2]);
    }
    simple_grads.push_back(grad);
  }

  assert( simple_grads.size() > 0 );
  for(unsigned int i = 0; i < simple_grads.size(); i++){
    assert( (simple_grads[i].size() + 1) == energy_bins.size() );
  }

  ice_gradient_array_pair_ = std::make_pair(energy_bins,simple_grads);
  DFWM.SetIceGradients(ice_gradient_array_pair_);

  if(!steeringParams_.quiet) std::cout << "Successfully loaded ice gradients." << std::endl;

  ice_gradient_loaded_ = true;
}

void GolemFit::LoadHoleIceResources(){
  holeIceSplines_.clear();
  if(!steeringParams_.quiet) std::cout << "Loading Hole Ice splines..." << std::endl;

  auto register_holeice_spline = [&](FluxComponent component, Topology topology) {
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
    holeIceSplines_.insert({std::make_pair(component,topology),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.holeice_spline_path+"/"+sample_name+"/"+sample_name+"_holeice_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits")))});
  };

  register_holeice_spline(FluxComponent::atmConv,Topology::track);
  register_holeice_spline(FluxComponent::atmPrompt,Topology::track);
  register_holeice_spline(FluxComponent::diffuseAstro_mu,Topology::track);

  if(steeringParams_.sampleToLoad == sampleTag::HESE){
    register_holeice_spline(FluxComponent::atmConv,Topology::shower);
    register_holeice_spline(FluxComponent::atmPrompt,Topology::shower);
    register_holeice_spline(FluxComponent::diffuseAstro_e,Topology::track);
    register_holeice_spline(FluxComponent::diffuseAstro_e,Topology::shower);
    register_holeice_spline(FluxComponent::diffuseAstro_mu,Topology::shower);
    register_holeice_spline(FluxComponent::diffuseAstro_tau,Topology::track);
    register_holeice_spline(FluxComponent::diffuseAstro_tau,Topology::shower);
  }

  assert(holeIceSplines_.size() != 0);

  // assuming all splines have the same validity ranges
  holeIceSpline_minValidValue_ = holeIceSplines_.begin()->second->lower_extent(2);
  holeIceSpline_maxValidValue_ = holeIceSplines_.begin()->second->upper_extent(2);

  if(!steeringParams_.quiet) std::cout << "Hole Ice splines forward parameter validity range set to [" + std::to_string(holeIceSpline_minValidValue_) + ", " + std::to_string(holeIceSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetHoleIceSplines(holeIceSplines_);

  holeice_resources_loaded_ = (true);
}

void GolemFit::LoadAttenuationResources(){
  attenuationSplines_.clear();
  if(!steeringParams_.quiet) std::cout << "Loading attenuation splines..." << std::endl;

  auto register_attenuation_spline = [&](FluxComponent component, LW::ParticleType particle_type) {
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
    attenuationSplines_.insert({{component, particle_type},
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.attenuation_spline_path+"/"+sample_name+"/"+sample_name+"_attenuation_spline_"+GetFluxComponentName(component)+"_"+GetParticleName(particle_type)+".fits")))});
  };

  register_attenuation_spline(FluxComponent::atmConv, LW::ParticleType::NuMu);
  register_attenuation_spline(FluxComponent::atmConv, LW::ParticleType::NuMuBar);
  register_attenuation_spline(FluxComponent::atmPrompt, LW::ParticleType::NuMu);
  register_attenuation_spline(FluxComponent::atmPrompt, LW::ParticleType::NuMuBar);
  register_attenuation_spline(FluxComponent::diffuseAstro_mu, LW::ParticleType::NuMu);
  register_attenuation_spline(FluxComponent::diffuseAstro_mu, LW::ParticleType::NuMuBar);

  assert(attenuationSplines_.size() != 0);

  // assuming all splines have the same validity ranges
  attenuationSpline_minValidValue_ = attenuationSplines_.begin()->second->lower_extent(2);
  attenuationSpline_maxValidValue_ = attenuationSplines_.begin()->second->upper_extent(2);

  if(!steeringParams_.quiet) std::cout << "Attenuation splines forward parameter validity range set to [" + std::to_string(attenuationSpline_minValidValue_) + ", " + std::to_string(attenuationSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetAttenuationSplines(attenuationSplines_);

  attenuation_resources_loaded_ = (true);
}

void GolemFit::LoadKaonEnergyLossesUncertaintiesResources(){
  if(!steeringParams_.quiet) std::cout << "Loading Kaon Energy Losses Uncertainty resources..." << std::endl;
  kaonLossesUncertaintySpline_ = std::make_shared<splinetable<>>(CheckedFilePath(dataPaths_.atmospheric_kaonlosses_spline_location));
  DFWM.SetKaonEnergyLossesUncertaintySpline(kaonLossesUncertaintySpline_);
  kaon_losses_spline_loaded_ = true;
}

void GolemFit::LoadAtmosphericDensityUncertaintiesResources(){
  if(!steeringParams_.quiet) std::cout << "Loading Atmospheric Density Uncertainty resources..."<< std::endl;
  atmosphericDensityUncertaintySpline_ = std::make_shared<splinetable<>>(CheckedFilePath(dataPaths_.atmospheric_density_spline_location+"/"+ steeringParams_.atmospheric_density_spline_filename));
  DFWM.SetAtmosphericDensityUncertaintySpline(atmosphericDensityUncertaintySpline_);
  atmospheric_density_spline_loaded_ = true;
}

void GolemFit::LoadBarrUncertaintiesResources(){
  barr_parameters_nusquids_gradients_.clear();
  std::string model_label = steeringParams_.sterile_model_label;
  //std::string model_label = std::to_string(new_physics_params_.index)+"_"+std::to_string(new_physics_params_.dm41sq)+"_"+std::to_string(new_physics_params_.th14)+"_"+std::to_string(new_physics_params_.th24)+"_"+std::to_string(new_physics_params_.th34)+"_"+std::to_string(new_physics_params_.del14)+"_"+std::to_string(new_physics_params_.del24);

  if(!steeringParams_.quiet) std::cout << "Loading Barr Uncertainty resources..."<< std::endl;

  auto register_barr_flux = [&](BarrParameter barr_param) {
    std::string sample_name = GetSampleName(steeringParams_.sampleToLoad);
    std::string barr_name = GetBarrParameterName(barr_param);
    barr_parameters_nusquids_gradients_.insert({barr_param,
        std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.barr_resources_location+"/"+ "barr_grad_" + barr_name + "_" + model_label + ".hdf5"))});
  };

  register_barr_flux(BarrParameter::WP);
  register_barr_flux(BarrParameter::WM);
  register_barr_flux(BarrParameter::ZP);
  register_barr_flux(BarrParameter::ZM);
  register_barr_flux(BarrParameter::YP);
  register_barr_flux(BarrParameter::YM);

  barr_resources_loaded_ = (true);
}

/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void GolemFit::ConstructCrossSectionWeighter(){
  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(CheckedFilePath(dataPaths_.diff_neutrino_cc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_antineutrino_cc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_neutrino_nc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_antineutrino_nc_xs_spline_path)));
  xs_weighter_constructed_=true;
}

void GolemFit::ConstructFluxWeighter(){
  if(steeringParams_.use_precalculated_nusquids_fluxes){
    if(steeringParams_.use_single_conventional_atmospheric_neutrino_file){
      fluxConv_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.conventional_nusquids_atmospheric_file), steeringParams_.use_nusquids_height_sampling);
    } else {
      fluxKaon_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.kaon_nusquids_atmospheric_file), steeringParams_.use_nusquids_height_sampling);
      fluxPion_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.pion_nusquids_atmospheric_file), steeringParams_.use_nusquids_height_sampling);
    }
    fluxPrompt_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.prompt_nusquids_atmospheric_file), steeringParams_.use_nusquids_height_sampling);
  } else {
    bool nugen_compatible = (steeringParams_.simType == MCType::NuGen);
    fluxPion_ = std::make_shared<LW::atmosNeutrinoFlux>(NewNuFlux::makeFlux(steeringParams_.atmospheric_modelName),nugen_compatible);
    std::dynamic_pointer_cast<NewNuFlux::PionKaonAdjustable>(std::dynamic_pointer_cast<LW::atmosNeutrinoFlux>(fluxPion_)->get())->setRelativeKaonContribution(0); //turn off kaons for the pion component
    std::dynamic_pointer_cast<NewNuFlux::KneeReweightable>(std::dynamic_pointer_cast<LW::atmosNeutrinoFlux>(fluxPion_)->get())->setKneeReweightingModel(steeringParams_.atmospheric_kneecorrection_modelName);

    fluxKaon_ = std::make_shared<LW::atmosNeutrinoFlux>(NewNuFlux::makeFlux(steeringParams_.atmospheric_modelName),nugen_compatible);
    std::dynamic_pointer_cast<NewNuFlux::PionKaonAdjustable>(std::dynamic_pointer_cast<LW::atmosNeutrinoFlux>(fluxKaon_)->get())->setRelativePionContribution(0); //turn off pion for kaon the component
    std::dynamic_pointer_cast<NewNuFlux::KneeReweightable>(std::dynamic_pointer_cast<LW::atmosNeutrinoFlux>(fluxKaon_)->get())->setKneeReweightingModel(steeringParams_.atmospheric_kneecorrection_modelName);

    fluxPrompt_ = std::make_shared<LW::atmosNeutrinoFlux>(NewNuFlux::makeFlux(steeringParams_.prompt_modelName),nugen_compatible);
    if(steeringParams_.correct_prompt_knee and steeringParams_.prompt_modelName == "sarcevic_std")
      std::dynamic_pointer_cast<NewNuFlux::KneeReweightable>(std::dynamic_pointer_cast<LW::atmosNeutrinoFlux>(fluxPrompt_)->get())->setKneeReweightingModel(steeringParams_.atmospheric_kneecorrection_modelName);
  }

  fluxAstro_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.astro_nusquids_file));
  //fluxAstro_ = std::make_shared<LW::PowerLawFlux>(steeringParams_.baseline_astro_normalization,
  //                                                steeringParams_.baseline_astro_spectral_index);

  fluxAstroGalactic_ = std::make_shared<LW::PowerLawFlux>(steeringParams_.baseline_galactic_astro_normalization,
                                                          steeringParams_.baseline_galactic_astro_spectral_index);

  flux_weighter_constructed_= true;
  simulation_initialized_ = false;
}

void GolemFit::ConstructMonteCarloGenerationWeighter(){
  mcw_.clear();
  std::map<std::string,MCSet> simInfo;
  if(steeringParams_.sampleToLoad == sampleTag::Sterile)
    simInfo = sterile::GetSimInfo(dataPaths_.mc_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + "/");
  else
    throw std::runtime_error("Simulation information not implemented for " + GetSampleName(steeringParams_.sampleToLoad) + "/");

  std::vector<MCSet> simSetsToLoad = sterile::GetSimulationSets(steeringParams_.simToLoad, dataPaths_.mc_path+ "/" + GetSampleName(steeringParams_.sampleToLoad));
  for(auto simset: simSetsToLoad ){
    for(auto g : simset.generators){
      mcw_.emplace_back(g);
    }
  }

  assert(not mcw_.empty());
  mc_generation_weighter_constructed_=true;
}

void GolemFit::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed_)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not xs_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  if(steeringParams_.use_single_conventional_atmospheric_neutrino_file){
    convFluxWeighter_ = std::make_shared<LW::Weighter>(fluxConv_,xsw_,mcw_);
  } else {
    pionFluxWeighter_ = std::make_shared<LW::Weighter>(fluxPion_,xsw_,mcw_);
    kaonFluxWeighter_ = std::make_shared<LW::Weighter>(fluxKaon_,xsw_,mcw_);
  }

  promptFluxWeighter_ = std::make_shared<LW::Weighter>(fluxPrompt_,xsw_,mcw_);
  astroFluxWeighter_ = std::make_shared<LW::Weighter>(fluxAstro_,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}

void GolemFit::ConstructOversizeWeighter(){
  if(steeringParams_.oversizeFunction == "NullCorrection")
    osw_ = OversizeWeighter();
  else
    osw_ = OversizeWeighter(CheckedFilePath(dataPaths_.oversize_path+"/"+steeringParams_.oversizeFunction+".dat"));
  oversize_weighter_constructed_=true;
}

/*************************************************************************************************************
 * Loading astrophysical templates
 * **********************************************************************************************************/

void GolemFit::LoadHealpixSkymap(){
  if(!steeringParams_.quiet) std::cout<<"Loading healpix skymap." <<std::endl;
  read_Healpix_map_from_fits<double>(CheckedFilePath(dataPaths_.skymaps_spline_path + "/" + steeringParams_.pi0_skymap),
                                     galactic_template_,1,2);
  int number_of_pixel = galactic_template_.Npix();
  double pixel_area = 4.*boost::math::constants::pi<double>()/number_of_pixel; // sr
  skymaps_galactic_template_normalization_ = 1./pixel_area;
  skymaps_loaded_=true;
}

void GolemFit::LoadSimpleHDF5Skymap(){
  if(!steeringParams_.quiet) std::cout<<"Loading simple format hdf5 skymap." <<std::endl;
  herr_t status;
  std::string str = CheckedFilePath(dataPaths_.skymaps_spline_path + "/" + steeringParams_.pi0_skymap);
  hid_t file_id = H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0)
    throw std::runtime_error("GolemFit::LoadSimpleHDF5Skymap::Error in trying to open file: " + str + ".");

  hid_t dataset = H5Dopen2(file_id, "flux", H5P_DEFAULT);

  // struct type in hdf5
  hid_t coso_id = H5Tcreate(H5T_COMPOUND, sizeof(double));
  status = H5Tinsert(coso_id,"value",0, H5T_NATIVE_DOUBLE);
  if(status<0)
    throw std::runtime_error("GolemFit::LoadSimpleHDF5Skymap::Error trying to construct type.");

  size_t size = H5Dget_storage_size(dataset)/sizeof(double);
  skymap_flux_array_.reset(new double[size]);
  status = H5Dread(dataset, coso_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, skymap_flux_array_.get());
  if(status<0)
    throw std::runtime_error("GolemFit::LoadSimpleHDF5Skymap::Error reading data set from file.");

  // close dataset
  H5Dclose(dataset);
  // close file
  H5Fclose(file_id);
  // close all HDF5 variables/memory
  H5close();
}

void GolemFit::LoadSkymaps(){
  LoadSimpleHDF5Skymap();
  skymaps_loaded_=true;
}

/*************************************************************************************************************
 * Functions to initialize the MC weights
 * **********************************************************************************************************/

void GolemFit::WeightMC(){
  if(not lepton_weighter_constructed_ and steeringParams_.simType==MCType::LeptonInjector)
    throw std::runtime_error("LeptonWeighter has to be constructed first.");
  if(not oversize_weighter_constructed_ and steeringParams_.simType==MCType::LeptonInjector)
    throw std::runtime_error("OversizeWeighter has to be constructed first.");
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");

  if(!simulation_initialized_)
    InitializeSimulationWeights();

  simulation_initialized_=true;
}

void GolemFit::InitializeSimulationWeights(){
  if(!steeringParams_.quiet)
    std::cout << "Initializing simulation weights" << std::endl;
  if(not dom_efficiency_splines_constructed_ and steeringParams_.spline_dom_efficiency)
    throw std::runtime_error("Simulation cannot be weighted until dom splines are loaded.");
  using iterator=std::deque<Event>::iterator;

  if(steeringParams_.simType == MCType::LeptonInjector){
    std::map<std::string,MCSet> simInfo;
    if(steeringParams_.sampleToLoad == sampleTag::Sterile)
      simInfo = sterile::GetSimInfo(dataPaths_.mc_path + "/" + GetSampleName(steeringParams_.sampleToLoad) + "/");
    else
      throw std::runtime_error("Simulation information not implemented for " + GetSampleName(steeringParams_.sampleToLoad) + "/");

    auto cache=[&](iterator it, iterator end){
      for(; it!=end; it++){
        auto& e=*it;
        LW::Event lw_e {
          e.primaryType,
          e.final_state_particle_0,
          e.final_state_particle_1,
          e.intX,
          e.intY,
          e.primaryEnergy,
          e.primaryAzimuth,
          e.primaryZenith,
          0.0,//e.x,
          0.0,//e.y,
          0.0,//e.z,
          0.0,//e.r,
          e.totalColumnDepth
        };

        //std::cout << e.intX << " " << e.intY << " " << e.primaryEnergy << " " << e.primaryAzimuth << " ";
        //std::cout << e.primaryZenith << e.x << " " << e.y << " " << e.z << " " << e.r << " " << e.totalColumnDepth << std::endl;

        double osweight = osw_.EvaluateOversizeCorrection(e.energy, e.zenith);
        try{
          if(steeringParams_.use_single_conventional_atmospheric_neutrino_file){
            e.oneWeight = convFluxWeighter_->get_oneweight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events; // you can use any "flux" weighter to get this one
          }else{
            e.oneWeight = pionFluxWeighter_->get_oneweight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events; // you can use any "flux" weighter to get this one
          }
          //e.cachedLivetime = 1.;
          //osweight = 1.;
          if(steeringParams_.use_single_conventional_atmospheric_neutrino_file){
            e.cachedConvWeight=(*convFluxWeighter_)(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            if(steeringParams_.add_contribution_from_tau_neutrinos){
              //std::cout << (*convFluxWeighter_)(lw_e) << " " << convFluxWeighter_->get_effective_tau_weight(lw_e) << " ";
              //std::cout << convFluxWeighter_->get_effective_tau_weight(lw_e)/(*convFluxWeighter_)(lw_e) << std::endl;
              e.cachedConvWeight+=convFluxWeighter_->get_effective_tau_weight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            }
            e.cachedConvKaonWeight = 0.0;
            e.cachedConvPionWeight = e.cachedConvWeight;
          } else {
            e.cachedConvPionWeight=(*pionFluxWeighter_)(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            if(steeringParams_.add_contribution_from_tau_neutrinos)
              e.cachedConvPionWeight+=pionFluxWeighter_->get_effective_tau_weight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            e.cachedConvKaonWeight=(*kaonFluxWeighter_)(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            if(steeringParams_.add_contribution_from_tau_neutrinos)
              e.cachedConvKaonWeight+=kaonFluxWeighter_->get_effective_tau_weight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
            e.cachedConvWeight = e.cachedConvPionWeight + e.cachedConvKaonWeight;
          }

          e.cachedPromptWeight=(*promptFluxWeighter_)(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
          if(steeringParams_.add_contribution_from_tau_neutrinos)
            e.cachedPromptWeight+=promptFluxWeighter_->get_effective_tau_weight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;

          if(e.primaryType == LW::ParticleType::NuMu or e.primaryType == LW::ParticleType::NuMuBar){
            e.cachedAstroMuWeight=(*astroFluxWeighter_)(lw_e)*(e.cachedLivetime*osweight/e.number_of_generated_mc_events);
            if(e.cachedAstroMuWeight == 0){
              std::cout << lw_e.energy << " " << lw_e.azimuth << " " << lw_e.zenith << std::endl;
              std::cout << (*fluxAstro_)(lw_e) << " " << (*xsw_)(lw_e) << " " << (*mcw_.front())(lw_e) << std::endl;
            }
            if(steeringParams_.add_contribution_from_tau_neutrinos)
            e.cachedAstroMuWeight+=astroFluxWeighter_->get_effective_tau_weight(lw_e)*e.cachedLivetime*osweight/e.number_of_generated_mc_events;
          }

          if(std::isnan(e.cachedConvPionWeight)){
            std::cout << "Atetmpting to weight event, but encountered nan. See:" << std::endl;
            std::cout << "Pion flux: " << (*pionFluxWeighter_)(lw_e) << std::endl;
            std::cout << "Lifetime: " << e.cachedLivetime*osweight << std::endl;
            std::cout << "Inverse of number of MC events: " << 1./e.number_of_generated_mc_events << std::endl;
            throw std::runtime_error("While attempting to initialize MC weights encounter nan.");
          }

          if(std::isnan(e.cachedConvKaonWeight)){
            std::cout << "Atetmpting to weight event, but encountered nan. See:" << std::endl;
            std::cout << "Kaon flux: " << (*kaonFluxWeighter_)(lw_e) << std::endl;
            std::cout << "Lifetime: " << e.cachedLivetime*osweight << std::endl;
            std::cout << "Inverse of number of MC events: " << 1./e.number_of_generated_mc_events << std::endl;
            throw std::runtime_error("While attempting to initialize MC weights encounter nan.");
          }

          if(std::isnan(e.cachedConvWeight)){
            std::cout << "Atetmpting to weight event, but encountered nan. See:" << std::endl;
            std::cout << "Conv flux: " << (*convFluxWeighter_)(lw_e) << std::endl;
            std::cout << "Lifetime: " << e.cachedLivetime*osweight << std::endl;
            std::cout << "Inverse of number of MC events: " << 1./e.number_of_generated_mc_events << std::endl;
            throw std::runtime_error("While attempting to initialize MC weights encounter nan.");
          }

          if(std::isnan(e.cachedPromptWeight)){
            std::cout << "Atetmpting to weight event, but encountered nan. See:" << std::endl;
            std::cout << "Prompt flux: " << (*promptFluxWeighter_)(lw_e) << std::endl;
            std::cout << "Lifetime: " << e.cachedLivetime*osweight << std::endl;
            std::cout << "Inverse of number of MC events: " << 1./e.number_of_generated_mc_events << std::endl;
            throw std::runtime_error("While attempting to initialize MC weights encounter nan.");
          }

          if(std::isnan(e.oneWeight)){
            std::cout << "Atetmpting to compute the one weight event, but encountered nan. See:" << std::endl;
            std::cout << e << std::endl;
            throw std::runtime_error("While attempting to initialize MC weights encounter nan.");
          }

          for(const auto& entry : barr_parameters_nusquids_gradients_){
            switch (entry.first) {
              case BarrParameter::WP:
                e.cachedBarrModWP = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::WM:
                e.cachedBarrModWM = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::YP:
                e.cachedBarrModYP = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::YM:
                e.cachedBarrModYM = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::ZP:
                e.cachedBarrModZP = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::ZM:
                e.cachedBarrModZM = (*(entry.second))(lw_e)*e.oneWeight;
                break;
              case BarrParameter::HP:
                // unused
                break;
              case BarrParameter::HM:
                // unused
                break;
              default:
                throw std::runtime_error("Impossible barr paramter label");
            }
          }
        } catch (std::exception exception){
          std::cout << exception.what() << " Encountered impossible to weight event. Removing it from deque. CAD." << std::endl;
          it = mainSimulation_.erase(it);
          /*
          e.oneWeight = 0;
          e.cachedConvPionWeight = 0;
          e.cachedConvKaonWeight = 0;
          e.cachedPromptWeight = 0;
          e.cachedBarrModWP = 0;
          e.cachedBarrModWM = 0;
          e.cachedBarrModYP = 0;
          e.cachedBarrModYM = 0;
          e.cachedBarrModZP = 0;
          e.cachedBarrModZM = 0;
          */
        }
      }
    };

    iterator it=mainSimulation_.begin(), end=mainSimulation_.end();
    cache(it,end);
  } else if (steeringParams_.simType == MCType::NuGen){
    throw std::runtime_error("Disable in ultra fast sterile mode: InitializeSimulationWeights");
  }
  else
    throw std::runtime_error("Invalid simulation type. No idea how to weight it. CA.");
}

void GolemFit::CalculateEarthTransferFunction(){
  if(steeringParams_.earth_propagation_interface == EarthPropagationInterfaces::nuSQuIDS)
    CalculateEarthTransferFunction_nuSQuIDS();
  else if (steeringParams_.earth_propagation_interface == EarthPropagationInterfaces::nuFATE)
    CalculateEarthTransferFunction_nuFATE();
  else
    throw std::runtime_error("GolemFit::Invalid Earth propagation interface.");
}

void GolemFit::CalculateEarthTransferFunction_nuFATE(){}

void GolemFit::CalculateEarthTransferFunction_nuSQuIDS(){}

double GolemFit::EvaluateEarthTransferFunction(Event & e, FluxComponent flux_component) const {
  if(steeringParams_.earth_propagation_interface == EarthPropagationInterfaces::nuSQuIDS) {
    //double energy = e.primaryEnergy*units.GeV;
    //double costh = cos(e.primaryZenith);
    //std::pair<unsigned int, unsigned int> nt = LW::Convert_PDG_Id_To_nuSQuIDS_Id(e.primaryType);
    switch(flux_component){
      // case FluxComponent::atmPion: return nus_atm_pion_.EvalFlavor(costh,energy,nt.first,nt.second);
      // case FluxComponent::atmKaon: return nus_atm_kaon_.EvalFlavor(costh,energy,nt.first,nt.second);
      // case FluxComponent::atmPrompt: return nus_atm_prompt_.EvalFlavor(costh,energy,nt.first,nt.second);
      // case FluxComponent::diffuseAstro: return nus_astro_.EvalFlavor(costh,energy,nt.first,nt.second);
      default: std::runtime_error("GolemFit::EvaluateEarthTransferFunction unimplemented component hit");
    }
  } else if (steeringParams_.earth_propagation_interface == EarthPropagationInterfaces::nuFATE)
    throw std::runtime_error("GolemFit::EvaluateEarthTransferFunction with nuFATE is not yet implemented.");
  else
    throw std::runtime_error("GolemFit::Invalid Earth propagation interface.");
  return std::numeric_limits<double>::quiet_NaN();
}

/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

void GolemFit::ConstructDataHistogram(){
  if(not data_loaded_)
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");

  // ALEJO TODO: maybe add separate boundaries for both histograms

  typedef std::remove_reference<decltype(std::get<0>(dataHist_))>::type Hist0;
  typedef std::remove_reference<decltype(std::get<0>(auxdataHist_))>::type Hist1;

  Hist0 h0(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                       LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth), // zenith dimension
                       LinearAxis(steeringParams_.azimuthbinEdge, steeringParams_.azimuthbinWidth), // azimuth dimension
                       LinearAxis(0,1), // topology dimension
                       LinearAxis(0,1)); // time dimension
  Hist1 h1(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                       LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth), // zenith dimension
                       LinearAxis(steeringParams_.azimuthbinEdge, steeringParams_.azimuthbinWidth), // azimuth dimension
                       LinearAxis(0,1), // topology dimension
                       LinearAxis(0,1)); // time dimension

  dataHist_ = std::make_tuple(h0);
  auxdataHist_ = std::make_tuple(h1);

  auto& data0 = std::get<0>(dataHist_);
  data0.getAxis(0)->setLowerLimit(steeringParams_.minFitEnergy);
  data0.getAxis(0)->setUpperLimit(steeringParams_.maxFitEnergy);
  data0.getAxis(1)->setLowerLimit(steeringParams_.minCosth);
  data0.getAxis(1)->setUpperLimit(steeringParams_.maxCosth);
  data0.getAxis(2)->setLowerLimit(steeringParams_.minAzimuth);
  data0.getAxis(2)->setUpperLimit(steeringParams_.maxAzimuth);

  auto& data1 = std::get<0>(auxdataHist_);
  data1.getAxis(0)->setLowerLimit(steeringParams_.minFitEnergy);
  data1.getAxis(0)->setUpperLimit(steeringParams_.maxFitEnergy);
  data1.getAxis(1)->setLowerLimit(steeringParams_.minCosth);
  data1.getAxis(1)->setUpperLimit(steeringParams_.maxCosth);
  data1.getAxis(2)->setLowerLimit(steeringParams_.minAzimuth);
  data1.getAxis(2)->setUpperLimit(steeringParams_.maxAzimuth);

  // fill in the histogram with the data
  bin(sample_, dataHist_, auxdataHist_, binner, steeringParams_.is_sideband);

  data_histogram_constructed_=true;
}

void GolemFit::ConstructSimulationHistogram(){
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before simulation histogram.");

  simHist_ = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(dataHist_)));
  auxsimHist_ = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(auxdataHist_)));
  bin(mainSimulation_, simHist_, auxsimHist_, binner, steeringParams_.is_sideband);

  simulation_histogram_constructed_=true;
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

hist_marray GolemFit::GetDataDistributionAll() const {
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before asking for it.");

  const auto& dataHist0 = std::get<0>(dataHist_);
  const auto& dataHist1 = std::get<0>(auxdataHist_);

  marray<double,5> array0 {static_cast<size_t>(dataHist0.getBinCount(4)),
                          static_cast<size_t>(dataHist0.getBinCount(3)),
                          static_cast<size_t>(dataHist0.getBinCount(2)),
                          static_cast<size_t>(dataHist0.getBinCount(1)),
                          static_cast<size_t>(dataHist0.getBinCount(0))};

  for(size_t iy=0; iy<dataHist0.getBinCount(4); iy++){ // year
    for(size_t it=0; it<dataHist0.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<dataHist0.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<dataHist0.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<dataHist0.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(dataHist0(ie,ic,ia,it,iy));
            array0[iy][it][ia][ic][ie] = 0;
            for(Event event : itc){
              array0[iy][it][ia][ic][ie] += event.cachedWeight;
            }
          }
        }
      }
    }
  }

  marray<double,5> array1 {static_cast<size_t>(dataHist1.getBinCount(4)),
                          static_cast<size_t>(dataHist1.getBinCount(3)),
                          static_cast<size_t>(dataHist1.getBinCount(2)),
                          static_cast<size_t>(dataHist1.getBinCount(1)),
                          static_cast<size_t>(dataHist1.getBinCount(0))};

  for(size_t iy=0; iy<dataHist1.getBinCount(4); iy++){ // year
    for(size_t it=0; it<dataHist1.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<dataHist1.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<dataHist1.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<dataHist1.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(dataHist1(ie,ic,ia,it,iy));
            array1[iy][it][ia][ic][ie] = 0;
            for(Event event : itc){
              array1[iy][it][ia][ic][ie] += event.cachedWeight;
            }
          }
        }
      }
    }
  }

  return std::make_tuple(array0, array1);
}

hist0_marray GolemFit::GetDataDistribution() const {return std::get<0>(GetDataDistributionAll());};

hist_marray GolemFit::GetRawExpectationAll(std::vector<double> fit_parameters) const {
  std::function<double(const Event &)> f = [](const Event & e){return 1.;};
  return GetWeightedExpectationAll(f);
}

hist0_marray GolemFit::GetRawExpectation(std::vector<double> fit_parameters) const {return std::get<0>(GetRawExpectationAll(fit_parameters));};

hist_marray GolemFit::GetExpectationAll(std::vector<double> fit_parameters) const {
  auto weighter = DFWM(fit_parameters);
  std::function<double(const Event&)> f = [&weighter](const Event & e){return weighter(e);};
  return GetWeightedExpectationAll(f);
}

hist0_marray GolemFit::GetExpectation(std::vector<double> fit_parameters) const {return std::get<0>(GetExpectationAll(fit_parameters));};

hist_marray GolemFit::GetSquareExpectationAll(std::vector<double> fit_parameters) const {
  auto weighter = DFWM(fit_parameters);
  std::function<double(const Event &)> f = [&weighter](const Event & e){return pow(weighter(e),2.);};
  return GetWeightedExpectationAll(f);
}

hist0_marray GolemFit::GetSquareExpectation(std::vector<double> fit_parameters) const {return std::get<0>(GetSquareExpectationAll(fit_parameters));};

hist_marray GolemFit::GetMaxExpectationErrorAll(FitParameters fit_parameters) const {
  return GetMaxExpectationErrorAll(ConvertFitParameters(fit_parameters));
}

hist0_marray GolemFit::GetMaxExpectationError(FitParameters fit_parameters) const {return std::get<0>(GetMaxExpectationErrorAll(fit_parameters));};

size_t GolemFit::GetTotalNumberOfBins() const {
  size_t bin_num_data = 1;
  const auto & dataHist1 = std::get<0>(dataHist_);
  for(unsigned int i = 0; i < dataHist1.getDimensions(); i++){
    bin_num_data *= dataHist1.getBinCount(i);
  }
  size_t bin_num_mc = 1;
  const auto & simHist1 = std::get<0>(simHist_);
  for(unsigned int i = 0; i < simHist1.getDimensions(); i++){
    bin_num_mc *= simHist1.getBinCount(i);
  }
  std::cout << bin_num_mc << "  "  << bin_num_data << std::endl;
 
  assert(bin_num_data == bin_num_mc);
  return bin_num_data;
}

namespace {
  struct square{
    template<typename T>
    void operator()(T& t) {
      for(auto & element : t)
        element = sqrt(element);
    }
};
}

hist_marray GolemFit::GetMaxExpectationErrorAll(std::vector<double> fit_parameters) const {
  auto sq_expectations = GetSquareExpectationAll(fit_parameters);
  square s;
  golemfit::fastmode::detail::apply_to_tuple(s, sq_expectations);
  return sq_expectations;
}
hist0_marray GolemFit::GetMaxExpectationError(std::vector<double> fit_parameters) const {return std::get<0>(GetMaxExpectationErrorAll(fit_parameters));};

nusquids::marray<std::vector<Event>,5> GolemFit::GetEventInSimulationHistogram() const {
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before asking for distributions.");
  const auto& simHist0 = std::get<0>(simHist_);

  marray<std::vector<Event>,5> event_array0 {static_cast<size_t>(simHist0.getBinCount(4)),
                                 static_cast<size_t>(simHist0.getBinCount(3)),
                                 static_cast<size_t>(simHist0.getBinCount(2)),
                                 static_cast<size_t>(simHist0.getBinCount(1)),
                                 static_cast<size_t>(simHist0.getBinCount(0))};

  for(size_t iy=0; iy<simHist0.getBinCount(4); iy++){ // year
    for(size_t it=0; it<simHist0.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<simHist0.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<simHist0.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<simHist0.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist0(ie,ic,ia,it,iy));
            //event_array0[iy][it][ia][ic][ie] = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist0(ie,ic,ia,it,iy)).entries();
            for(auto event : itc.entries()){
             event_array0[iy][it][ia][ic][ie].push_back(event);
            }
          }
        }
      }
    }
  }

  return event_array0;
}

hist_marray GolemFit::GetWeightedExpectationAll(std::function<double(const Event &)> f) const {
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before asking for distributions.");

  const auto& simHist0 = std::get<0>(simHist_);
  const auto& simHist1 = std::get<0>(auxsimHist_);

  marray<double,5> array0 {static_cast<size_t>(simHist0.getBinCount(4)),
                          static_cast<size_t>(simHist0.getBinCount(3)),
                          static_cast<size_t>(simHist0.getBinCount(2)),
                          static_cast<size_t>(simHist0.getBinCount(1)),
                          static_cast<size_t>(simHist0.getBinCount(0))};
  std::fill(array0.begin(),array0.end(),0);

  for(size_t iy=0; iy<simHist0.getBinCount(4); iy++){ // year
    for(size_t it=0; it<simHist0.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<simHist0.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<simHist0.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<simHist0.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist0(ie,ic,ia,it,iy));
            double expectation=0;
            for(auto event : itc.entries()){
              expectation+=f(event);
            }
            assert(expectation>=0.0 && "Expectation cannot be negative");

            array0[iy][it][ia][ic][ie] = expectation;
          }
        }
      }
    }
  }
 
  marray<double,5> array1 {static_cast<size_t>(simHist1.getBinCount(4)),
                          static_cast<size_t>(simHist1.getBinCount(3)),
                          static_cast<size_t>(simHist1.getBinCount(2)),
                          static_cast<size_t>(simHist1.getBinCount(1)),
                          static_cast<size_t>(simHist1.getBinCount(0))};
  std::fill(array1.begin(),array1.end(),0);

  for(size_t iy=0; iy<simHist1.getBinCount(4); iy++){ // year
    for(size_t it=0; it<simHist1.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<simHist1.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<simHist1.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<simHist1.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist1(ie,ic,ia,it,iy));
            double expectation=0;
            for(auto event : itc.entries()){
              expectation+=f(event);
            }
            assert(expectation>=0.0 && "Expectation cannot be negative");

            array1[iy][it][ia][ic][ie] = expectation;
          }
        }
      }
    }
  }

  return std::make_tuple(array0, array1);
}

hist0_marray GolemFit::GetWeightedExpectation(std::function<double(const Event &)> f) const {return std::get<0>(GetWeightedExpectationAll(f));};

hist_marray GolemFit::GetExpectationAll(FitParameters fit_params) const {
  return GetExpectationAll(ConvertFitParameters(fit_params));
}
hist0_marray GolemFit::GetExpectation(FitParameters fit_params) const {return std::get<0>(GetExpectationAll(fit_params));};

hist_marray GolemFit::GetSquareExpectationAll(FitParameters fit_params) const {
  return GetSquareExpectationAll(ConvertFitParameters(fit_params));
}
hist0_marray GolemFit::GetSquareExpectation(FitParameters fit_params) const {return std::get<0>(GetSquareExpectationAll(fit_params));};

hist_marray GolemFit::GetRealizationAll(std::vector<double> fit_params, int seed) const {
  std::mt19937 rng;
  rng.seed(seed);

  if(!steeringParams_.quiet)
    std::cout << "construct weighter" << std::endl;
  auto weighter=DFWM(fit_params);

  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    if(std::isnan(w) || std::isinf(w) || w<0){
      std::cout << "Bad weight! " << w << std::endl;
      std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
      std::cout << e.energy << ' ' << cos(e.zenith) <<' ' << e.year << ' ' << w << std::endl;
    }
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization=phys_tools::likelihood::generateSample(weights,mainSimulation_,expected,rng);
  auto fullRealizationHist = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(dataHist_)));
  auto auxfullRealizationHist = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(auxdataHist_)));
  bin(realization,fullRealizationHist,auxfullRealizationHist,binner, steeringParams_.is_sideband);
  auto& realizationHist0 = std::get<0>(fullRealizationHist);
  auto& realizationHist1 = std::get<0>(auxfullRealizationHist);

  if(realization.size() == 0){
    throw std::runtime_error("No events generated. Expected events are "+std::to_string(expected));
  }

  marray<double,5> array0 {static_cast<size_t>(realizationHist0.getBinCount(4)),
                          static_cast<size_t>(realizationHist0.getBinCount(3)),
                          static_cast<size_t>(realizationHist0.getBinCount(2)),
                          static_cast<size_t>(realizationHist0.getBinCount(1)),
                          static_cast<size_t>(realizationHist0.getBinCount(0))};
  std::fill(array0.begin(),array0.end(),0);

  for(size_t iy=0; iy<realizationHist0.getBinCount(4); iy++){ // year
    for(size_t it=0; it<realizationHist0.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<realizationHist0.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<realizationHist0.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<realizationHist0.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(realizationHist0(ie,ic,ia,it,iy));
            array0[iy][it][ia][ic][ie] = itc.size();
          }
        }
      }
    }
  }

  marray<double,5> array1 {static_cast<size_t>(realizationHist1.getBinCount(4)),
                          static_cast<size_t>(realizationHist1.getBinCount(3)),
                          static_cast<size_t>(realizationHist1.getBinCount(2)),
                          static_cast<size_t>(realizationHist1.getBinCount(1)),
                          static_cast<size_t>(realizationHist1.getBinCount(0))};
  std::fill(array1.begin(),array1.end(),0);

  for(size_t iy=0; iy<realizationHist1.getBinCount(4); iy++){ // year
    for(size_t it=0; it<realizationHist1.getBinCount(3); it++){ // topology
      for(size_t ia=0; ia<realizationHist1.getBinCount(2); ia++){ // azimuth
        for(size_t ic=0; ic<realizationHist1.getBinCount(1); ic++){ // zenith
          for(size_t ie=0; ie<realizationHist1.getBinCount(0); ie++){ // energy
            auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(realizationHist1(ie,ic,ia,it,iy));
            array1[iy][it][ia][ic][ie] = itc.size();
          }
        }
      }
    }
  }

  return std::make_tuple(array0, array1);
}
hist0_marray GolemFit::GetRealization(std::vector<double> fit_params, int seed) const {return std::get<0>(GetRealizationAll(fit_params, seed));};

hist_marray GolemFit::GetRealizationAll(FitParameters nuisance, int seed) const {
  return GetRealizationAll(ConvertFitParameters(nuisance),seed);
}
hist0_marray GolemFit::GetRealization(FitParameters nuisance, int seed) const {return std::get<0>(GetRealizationAll(nuisance, seed));};

histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GolemFit::GetEnergyZenithExpectationHistogram(std::vector<double> fit_params, int topology) const {
  auto EnergyZenithHistogram = histogram<2,entryStoringBin<std::reference_wrapper<const Event>>>(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                                                                                                 LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth)); // zenith dimension

  for(auto & e : mainSimulation_){
    if(topology<0)
      EnergyZenithHistogram.add(e.energy,cos(e.zenith),amount(std::cref(e)));
    else if (e.topology == (unsigned int)(topology))
      EnergyZenithHistogram.add(e.energy,cos(e.zenith),amount(std::cref(e)));
  }

  return EnergyZenithHistogram;
}

histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GolemFit::GetEnergyZenithExpectationHistogram(FitParameters fit_params, int topology) const {
  return GetEnergyZenithExpectationHistogram(ConvertFitParameters(fit_params),topology);
}

std::function<double(const Event&)> GolemFit::GetEventWeighter(FitParameters fp) const {
  return DFWM(ConvertFitParameters(fp));
}

/*************************************************************************************************************
 * Functions to construct likelihood problem and evaluate it
 * **********************************************************************************************************/

void GolemFit::ConstructLikelihoodProblem(Priors pr, std::vector<FitParameters> nuisanceSeed, FitParametersFlag fixedParams){
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before likelihood problem can be formulated.");
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before likelihood problem can be formulated.");
  if(steeringParams_.spline_dom_efficiency and not dom_efficiency_splines_constructed_)
    throw std::runtime_error("DOM efficiency is going to be splinned, but you have not loaded any valid splines.");
  if(steeringParams_.spline_hole_ice and not holeice_resources_loaded_)
    throw std::runtime_error("Hole ice is going to be splinned, but you have not loaded any valid splines.");

  fixedParams_=fixedParams;
  fitSeed_=nuisanceSeed;
  ForceFitSeedSanity();

  // standard flux parameters
  LimitedGaussianPrior convNormPrior(pr.convNormCenter,pr.convNormWidth,1.e-5,std::numeric_limits<double>::max());
  LimitedGaussianPrior promptNormPrior(pr.promptNormCenter,pr.promptNormWidth,0.0,std::numeric_limits<double>::max());
  GaussianPrior crSlopePrior(pr.crSlopeCenter,pr.crSlopeWidth);
  GaussianPrior kaonPrior(pr.piKRatioCenter,pr.piKRatioWidth);
  GaussianPrior nanPrior(pr.nuNubarRatioCenter,pr.nuNubarRatioWidth);
  //GaussianPrior ZCPrior(0.0,pr.zenithCorrectionMultiplier*GetZenithCorrectionScale());
  GaussianPrior ZCPrior(pr.zenithCorrectionCenter,pr.zenithCorrectionWidth);

  // standard flux parameters
  LimitedGaussianPrior BarrWPPrior(pr.barrWPCenter,pr.barrWPWidth,pr.barrWPMin,pr.barrWPMax);
  LimitedGaussianPrior BarrWMPrior(pr.barrWMCenter,pr.barrWMWidth,pr.barrWMMin,pr.barrWMMax);
  LimitedGaussianPrior BarrYPPrior(pr.barrYPCenter,pr.barrYPWidth,pr.barrYPMin,pr.barrYPMax);
  LimitedGaussianPrior BarrYMPrior(pr.barrYMCenter,pr.barrYMWidth,pr.barrYMMin,pr.barrYMMax);
  LimitedGaussianPrior BarrZPPrior(pr.barrZPCenter,pr.barrZPWidth,pr.barrZPMin,pr.barrZPMax);
  LimitedGaussianPrior BarrZMPrior(pr.barrZMCenter,pr.barrZMWidth,pr.barrZMMin,pr.barrZMMax);
  LimitedGaussianPrior BarrHPPrior(pr.barrHPCenter,pr.barrHPWidth,pr.barrHPMin,pr.barrHPMax);
  LimitedGaussianPrior BarrHMPrior(pr.barrHMCenter,pr.barrHMWidth,pr.barrHMMin,pr.barrHMMax);

  // detector systematics
  LimitedGaussianPrior domEffPrior(pr.domEffCenter,pr.domEffWidth,domEffSpline_minValidValue_,domEffSpline_maxValidValue_);
  LimitedGaussianPrior hqdomEffPrior(pr.hqdomEffCenter,pr.hqdomEffWidth,hqdomEffSpline_minValidValue_,hqdomEffSpline_maxValidValue_);
  LimitedGaussianPrior holeiceForwardPrior(pr.holeiceForwardCenter,pr.holeiceForwardWidth,holeIceSpline_minValidValue_,holeIceSpline_maxValidValue_);
  LimitedGaussianPrior icegrad0Prior(pr.icegrad0Center,std::numeric_limits<double>::max(),pr.icegrad0Min,pr.icegrad0Max);
  LimitedGaussianPrior icegrad1Prior(pr.icegrad1Center,std::numeric_limits<double>::max(),pr.icegrad1Min,pr.icegrad1Max);
  LimitedGaussianPrior icegrad2Prior(pr.icegrad2Center,pr.icegrad2Width,pr.icegrad2Min,pr.icegrad2Max);

  LimitedGaussianPrior astroNormPrior(pr.astroNormCenter,pr.astroNormWidth,0.0,std::numeric_limits<double>::max());
  GaussianPrior astroDeltaGammaPrior(pr.astroDeltaGammaCenter,pr.astroDeltaGammaWidth);

  LimitedGaussianPrior astroNormSecPrior(pr.astroNormSecCenter,pr.astroNormSecWidth,0.0,std::numeric_limits<double>::max());
  GaussianPrior astroDeltaGammaSecPrior(pr.astroDeltaGammaSecCenter,pr.astroDeltaGammaSecWidth);

  GaussianPrior nuCrossPrior(pr.nuXSCenter,pr.nuXSWidth);
  GaussianPrior nuBarCrossPrior(pr.nubarXSCenter,pr.nubarXSWidth);

  GaussianPrior kaonLossesPrior(pr.kaonLossesCenter,pr.kaonLossesWidth);

  auto basicpriors=makePriorSet(
            convNormPrior, // conventional normalization prior
			      promptNormPrior, // prompt normalization prior
			      crSlopePrior, // CR ray slope prior
			      kaonPrior, // pik ratio prior
            nanPrior, // conventional particle balance
			      ZCPrior, // jordi-delta prior
			      BarrWPPrior, // barr parameter prior
			      BarrWMPrior, // barr parameter prior
			      BarrYPPrior, // barr parameter prior
			      BarrYMPrior, // barr parameter prior
			      BarrZPPrior, // barr parameter prior
			      BarrZMPrior, // barr parameter prior
			      BarrHPPrior, // barr parameter prior
			      BarrHMPrior, // barr parameter prior
			      domEffPrior, // dom efficiency prior
			      hqdomEffPrior, // dom efficiency prior
			      holeiceForwardPrior,// hole ice forward prior
			      icegrad0Prior, // ice gradient prior
			      icegrad1Prior, // ice gradient prior
			      icegrad2Prior, // ice gradient prior
            astroNormPrior, // astro norm prior
            astroDeltaGammaPrior, // astro delta gamma prior
            nuCrossPrior, // neutrino xs prior
            nuBarCrossPrior, // antineutrino xs prior
            kaonLossesPrior, // kaon losses prior
            astroNormSecPrior, // second astrophysical component normalization prior
            astroDeltaGammaSecPrior // second astrophysical component parameter prior
  );

  Gaussian2DPrior ice_gradient_joined_prior(pr.icegrad0Center,pr.icegrad1Center, pr.icegrad0Width, pr.icegrad1Width, pr.icegrad01_correlation);

  Gaussian2DPrior astro_correlated_prior(pr.astro1Comp2DNormCenter, pr.astro1Comp2DDeltaGammaCenter, pr.astro1Comp2DNormWidth, pr.astro1Comp2DDeltaGammaWidth, pr.astro1Comp2DCorrelation);

  auto llhpriors = makeArbitraryPriorSet<PriorIndices>(basicpriors,ice_gradient_joined_prior,astro_correlated_prior);

  auto fitseedvec = ConvertFitParameters(nuisanceSeed.front());

  prob_ = std::make_shared<LType>(phys_tools::likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>,27>(
                                  dataHist_, {simHist_}, llhpriors, {0.0}, simpleLocalDataWeighterConstructor(), DFWM,
                                  //phys_tools::likelihood::poissonLikelihood(), fitseedvec ));
                                  phys_tools::likelihood::SAYLikelihood(), fitseedvec));
  prob_->setEvaluationThreadCount(steeringParams_.evalThreads);

  auxprob_ = std::make_shared<LType>(phys_tools::likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>,27>(
                                  auxdataHist_, {auxsimHist_}, llhpriors, {0.0}, simpleLocalDataWeighterConstructor(), DFWM,
                                  //phys_tools::likelihood::poissonLikelihood(), fitseedvec ));
                                  phys_tools::likelihood::SAYLikelihood(), fitseedvec));
  auxprob_->setEvaluationThreadCount(steeringParams_.evalThreads);

  likelihood_problem_constructed_=true;
}

// look up zenith correction scale for a particular flux model.
double GolemFit::GetZenithCorrectionScale() const {
  std::map<std::string,double> delta_alpha {
    {"honda2006",8./7.},
    {"CombinedGHandHG_H3a_QGSJET",4. /7.},
    {"CombinedGHandHG_H3a_SIBYLL2",8./7.},
    {"PolyGonato_QGSJET-II-04",0.5},
    {"PolyGonato_SIBYLL2",1.0},
    {"ZatsepinSokolskaya_pamela_QGSJET",5./7.},
    {"ZatsepinSokolskaya_pamela_SIBYLL2",5./7.},
                  };
  if( delta_alpha.find(steeringParams_.atmospheric_modelName) == delta_alpha.end() )
    throw std::runtime_error("Jordi delta key not found. Aborting.");
  return delta_alpha[steeringParams_.atmospheric_modelName];
}

double GolemFit::EvalLLH(std::vector<double> nuisance, bool include_prior) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  return -prob_->evaluateLikelihood(nuisance,include_prior);
}

double GolemFit::EvalLLH(FitParameters nuisance, bool include_prior) const {
  return EvalLLH(ConvertFitParameters(nuisance),include_prior);
}

phys_tools::autodiff::FD<27> GolemFit::EvalLLHGradient(std::vector<phys_tools::autodiff::FD<27>> v) const {
  return -prob_->evaluateLikelihood(v);
}

void GolemFit::ForceFitSeedSanity(){
  for(auto & fitSeed: fitSeed_)
    ForceFitSeedSanity(fitSeed);
}

void GolemFit::ForceFitSeedSanity(FitParameters& fitSeed) {
  if(steeringParams_.spline_hole_ice and (fitSeed.holeiceForward < holeIceSpline_minValidValue_ or fitSeed.holeiceForward > holeIceSpline_maxValidValue_)){
    double proposed_holeice_seed = (holeIceSpline_maxValidValue_ + holeIceSpline_minValidValue_)/2.;
    if(!steeringParams_.quiet) std::cout << "Holeice forward seed value ("+std::to_string(fitSeed.holeiceForward)+ ") " + " outside of spline valid range [" + std::to_string(holeIceSpline_minValidValue_) + ", " + std::to_string(holeIceSpline_maxValidValue_) + "]. Changing it to " + std::to_string(proposed_holeice_seed) + "." << std::endl;
    fitSeed.holeiceForward = proposed_holeice_seed;
  }

  if(steeringParams_.spline_dom_efficiency and (fitSeed.domEfficiency < domEffSpline_minValidValue_ or fitSeed.domEfficiency > domEffSpline_maxValidValue_)){
    double proposed_dom_efficiency_seed = (domEffSpline_maxValidValue_ + domEffSpline_minValidValue_)/2.;
    if(!steeringParams_.quiet) std::cout << "DOM efficiency seed value ("+std::to_string(fitSeed.domEfficiency)+ ") " + " outside of spline valid range [" + std::to_string(domEffSpline_minValidValue_) + ", " + std::to_string(domEffSpline_maxValidValue_) + "]. Changing it to " + std::to_string(proposed_dom_efficiency_seed) + "." << std::endl;
    fitSeed.domEfficiency = proposed_dom_efficiency_seed;
  }

  if(steeringParams_.spline_hqdom_efficiency and (fitSeed.hqdomEfficiency < hqdomEffSpline_minValidValue_ or fitSeed.hqdomEfficiency > hqdomEffSpline_maxValidValue_)){
    double proposed_dom_efficiency_seed = (hqdomEffSpline_maxValidValue_ + hqdomEffSpline_minValidValue_)/2.;
    if(!steeringParams_.quiet) std::cout << "HQ DOM efficiency seed value ("+std::to_string(fitSeed.hqdomEfficiency)+ ") " + " outside of spline valid range [" + std::to_string(hqdomEffSpline_minValidValue_) + ", " + std::to_string(hqdomEffSpline_maxValidValue_) + "]. Changing it to " + std::to_string(proposed_dom_efficiency_seed) + "." << std::endl;
    fitSeed.hqdomEfficiency = proposed_dom_efficiency_seed;
  }
}

FitResult GolemFit::MinLLH(bool evaluate_sideband) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");

  if(!CheckAsimovSanity()) {
      throw std::runtime_error("Asimov set up improperly");
  }

  FitResult final_result;
  final_result.likelihood = std::numeric_limits<double>::max();

  for(auto fitSeed : fitSeed_){
    std::vector<double> seed=ConvertFitParameters(fitSeed);
    prob_->setSeed(seed);

    std::vector<unsigned int> fixedIndices;
    std::vector<bool> FixVec=ConvertFitParametersFlag(fixedParams_);
    for(size_t i=0; i!=FixVec.size(); i++)
        if(FixVec[i]) fixedIndices.push_back(i);

    phys_tools::lbfgsb::LBFGSB_Driver minimizer;
    minimizer.setGradientTolerance(1e-20);
    minimizer.setChangeTolerance(1e-20);

    //std::cout << domEffSpline_minValidValue_ << " " << domEffSpline_maxValidValue_ << std::endl;
    //std::cout << holeIceSpline_minValidValue_ << " " << holeIceSpline_maxValidValue_ << std::endl;

    // parameter seed, gradient factor, min boundary, max boundary
    minimizer.addParameter(seed[0],.001,1.e-5,10); // conv norm
    minimizer.addParameter(seed[1],.0001,0.0); // prompt norm
    minimizer.addParameter(seed[2],.001,-0.8,0.8); // cr slope
    minimizer.addParameter(seed[3],.001,0.0,2.0); // pik
    minimizer.addParameter(seed[4],.001, 0.0, 2.0); // conv particle balance
    minimizer.addParameter(seed[5],.001,-10.0,10.0); // jordi delta
    minimizer.addParameter(seed[6],.001,priors_.barrWPMin*0.95,priors_.barrWPMax*0.95); // barr wp
    minimizer.addParameter(seed[7],.001,priors_.barrWMMin*0.95,priors_.barrWMMax*0.95); // barr wm
    minimizer.addParameter(seed[8],.001,priors_.barrYPMin*0.95,priors_.barrYPMax*0.95); // barr yp
    minimizer.addParameter(seed[9],.001,priors_.barrYMMin*0.95,priors_.barrYMMax*0.95); // barr ym
    minimizer.addParameter(seed[10],.001,priors_.barrZPMin*0.95,priors_.barrZPMax*0.95); // barr zp
    minimizer.addParameter(seed[11],.001,priors_.barrZMMin*0.95,priors_.barrZMMax*0.95); // barr zm
    minimizer.addParameter(seed[12],.001,priors_.barrHPMin*0.95,priors_.barrHPMax*0.95); // barr hp
    minimizer.addParameter(seed[13],.001,priors_.barrHMMin*0.95,priors_.barrHMMax*0.95); // barr hm
    minimizer.addParameter(seed[14],.005,domEffSpline_minValidValue_ + 0.1*(domEffSpline_maxValidValue_ - domEffSpline_minValidValue_),
                                         domEffSpline_maxValidValue_ - 0.1*(domEffSpline_maxValidValue_ - domEffSpline_minValidValue_)); // dom eff
    minimizer.addParameter(seed[15],.005,hqdomEffSpline_minValidValue_ + 0.1*(hqdomEffSpline_maxValidValue_ - hqdomEffSpline_minValidValue_),
                                         hqdomEffSpline_maxValidValue_ - 0.1*(hqdomEffSpline_maxValidValue_ - hqdomEffSpline_minValidValue_)); // hq dom eff
    minimizer.addParameter(seed[16],.001,holeIceSpline_minValidValue_ + 0.1*(holeIceSpline_maxValidValue_ - holeIceSpline_minValidValue_),
                                         holeIceSpline_maxValidValue_ - 0.1*(holeIceSpline_maxValidValue_ - holeIceSpline_minValidValue_));// hole ice forward
    minimizer.addParameter(seed[17],.0005,-4,4); // ice grad 0
    minimizer.addParameter(seed[18],.0005,-4,4); // ice grad 1
    minimizer.addParameter(seed[19],.0005,-4,4); // ice grad 2
    minimizer.addParameter(seed[20],.0005,0.0); // astro norm
    minimizer.addParameter(seed[21],.0001,-5,5); // astro delta gamma
    minimizer.addParameter(seed[22],.001,attenuationSpline_minValidValue_ + 0.1*(attenuationSpline_maxValidValue_ - attenuationSpline_minValidValue_),
                                         attenuationSpline_maxValidValue_ - 0.1*(attenuationSpline_maxValidValue_ - attenuationSpline_minValidValue_));   // nuxs
    minimizer.addParameter(seed[23],.001,attenuationSpline_minValidValue_ + 0.1*(attenuationSpline_maxValidValue_ - attenuationSpline_minValidValue_),
                                         attenuationSpline_maxValidValue_ - 0.1*(attenuationSpline_maxValidValue_ - attenuationSpline_minValidValue_));   // nubarxs
    minimizer.addParameter(seed[24],.001,-1.0,1.0); // kaon losses
    minimizer.addParameter(seed[25],.001,0.0); // second astro component norm
    minimizer.addParameter(seed[26],.001,-2.,2.);// second astro component parameters

    minimizer.setHistorySize(20);

    //std::cout << "fix parameters" << std::endl;
    for(auto idx : fixedIndices){
      minimizer.fixParameter(idx);
      //std::cout << idx << std::endl;
    }

    FitResult result;
    result.succeeded=DoFitLBFGSB(*prob_, minimizer);
    result.likelihood=minimizer.minimumValue();
    result.params=ConvertVecToFitParameters(minimizer.minimumPosition());

    /*
    for(auto p : minimizer.minimumPosition()){
      std::cout << p << std::endl;
    }
    */

    result.nEval+=minimizer.numberOfEvaluations();
    result.nGrad+=minimizer.numberOfEvaluations();

    //std::cout << "Minimizer error message: " << minimizer.errorMessage() << std::endl;
    //std::cout << "Likelihood of baby fit : " << result.likelihood << std::endl;
    if(result.likelihood < final_result.likelihood)
      final_result = result;

    if(evaluate_sideband)
      final_result.aux_likelihood = -auxprob_->evaluateLikelihood(ConvertFitParameters(result.params),false);
  }

  return final_result;
}

/*************************************************************************************************************
 * Functions to apply new physics hypothesis
 * **********************************************************************************************************/

void GolemFit::SetNewPhysicsParams(NewPhysicsParams new_physics_params){
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot weight to sterile hypothesis without simulation.");

  new_physics_params_=new_physics_params;
  if(steeringParams_.simType == MCType::LeptonInjector){
    if(not mc_generation_weighter_constructed_)
      throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
    if(not xs_weighter_constructed_)
      throw std::runtime_error("Cross section weighter has to be constructed first.");
    ConstructFluxWeighter();
    ConstructLeptonWeighter();
  }
  WeightMC();
}

/*************************************************************************************************************
 * Implementation of the fast mode
 * **********************************************************************************************************/

void GolemFit::ConstructFastMode() {
  if((not simulation_loaded_ ) or mainSimulation_.empty())
    throw std::runtime_error("No simulation has been loaded. Cannot caress what does not exist.");
  if(not data_loaded_)
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");
  // This piece of code is inspired in ideas implemented in LVTools by CA. A better implementation
  // was done, independently, by BJPJ on the Sterilizer.

  double meta_scaling = steeringParams_.fastmode_scaling;
  if(not steeringParams_.quiet)
    std::cout << "Using metascaling " << meta_scaling << " large metascaling can produce inaccurate results. Love CAD." << std::endl;
  if((meta_scaling > 1.) and steeringParams_.paranoid)
    throw std::runtime_error("Paranoid mode is on. Using metascaling greater than one can give bad results. Reconsider.");

  using MetaHistType = histogram<4,entryStoringBin<std::reference_wrapper<const Event>>>;

  MetaHistType metaHist(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth*meta_scaling), // energy dimension
                        LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth*meta_scaling), // zenith dimension
                        LinearAxis(steeringParams_.azimuthbinEdge, steeringParams_.azimuthbinWidth), // azimuth dimension
                        LinearAxis(0,1)); // particle type dimension

  auto meta_binner = [](MetaHistType& h, const Event& e){
    h.add(e.primaryEnergy,cos(e.primaryZenith),e.primaryAzimuth,int(e.primaryType),amount(std::cref(e)));
  };

  struct combiner {
      Event operator()(const std::vector<std::reference_wrapper<const Event>>& events) {
        return combine_events(events);
      }
  };

  metaEvents_ = golemfit::fastmode::get_fastmode_events<Event>(metaHist, meta_binner, combiner(), simHist_);

  if(!steeringParams_.quiet)
    std::cout << "Went from this many MC events " << mainSimulation_.size() << " to " << metaEvents_.size() << std::endl;
  assert(!metaEvents_.empty());

  // reset meta events detector systematic properties
  if(steeringParams_.fastmode_detector_systematics_reset){
    if(steeringParams_.spline_dom_efficiency){
      if(domefficiencySplines_.empty())
        throw std::runtime_error("ConstructFastMode: DOM efficiency splines are empty. Cannot set cached efficiency.");
      for(auto & e : metaEvents_){
        domEffSetter_->setCache(e);
      }
    }
    if(steeringParams_.spline_hole_ice){
      if(holeIceSplines_.empty())
        throw std::runtime_error("ConstructFastMode: Hole ice splines are empty. Cannot set cached efficiency.");
      for(auto & e : metaEvents_){
        holeIceSetter_->setCache(e);
      }
    }
  }

  simHist_ = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(dataHist_)));
  auxsimHist_ = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(auxdataHist_)));
  bin(metaEvents_, simHist_, auxsimHist_, binner, steeringParams_.is_sideband);

  fastmode_constructed_ = true;
}

/*************************************************************************************************************
 * Functions to set options in the class
 * **********************************************************************************************************/

// Check that the directories where files are mean to be exist
bool GolemFit::CheckDataPaths(DataPaths dp) const {
  if(steeringParams_.simType == MCType::LeptonInjector){
    CheckDataPath(dp.squids_files_path);
    CheckDataPath(dp.prompt_squids_files_path);
    CheckedFilePath(dp.neutrino_cc_xs_spline_path);
    CheckedFilePath(dp.antineutrino_cc_xs_spline_path);
    CheckedFilePath(dp.neutrino_nc_xs_spline_path);
    CheckedFilePath(dp.antineutrino_nc_xs_spline_path);
    CheckDataPath(dp.oversize_path);
  }
  if(steeringParams_.spline_dom_efficiency)
    CheckDataPath(dp.domeff_spline_path);
  CheckDataPath(dp.compact_file_path);
  CheckDataPath(dp.data_path);
  CheckDataPath(dp.mc_path);
  return true;
}

// Check a directory exists and throw a relevant error otherwise.
bool GolemFit::CheckDataPath(std::string p) const {
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
   }
 else{
   std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
   return false;
 }
 return status;
}

// Given a human readable nuisance parameter set, make a nuisance vector
std::vector<double> GolemFit::ConvertFitParameters(FitParameters ns) const {
  std::vector<double> nuis;

  nuis.push_back(ns.convNorm);
  nuis.push_back(ns.promptNorm);
  nuis.push_back(ns.CRDeltaGamma);
  nuis.push_back(ns.piKRatio);
  nuis.push_back(ns.NeutrinoAntineutrinoRatio);
  nuis.push_back(ns.zenithCorrection);

  nuis.push_back(ns.barrWP);
  nuis.push_back(ns.barrWM);
  nuis.push_back(ns.barrYP);
  nuis.push_back(ns.barrYM);
  nuis.push_back(ns.barrZP);
  nuis.push_back(ns.barrZM);
  nuis.push_back(ns.barrHP);
  nuis.push_back(ns.barrHM);

  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.hqdomEfficiency);
  nuis.push_back(ns.holeiceForward);
  nuis.push_back(ns.icegrad0);
  nuis.push_back(ns.icegrad1);
  nuis.push_back(ns.icegrad2);
  nuis.push_back(ns.astroNorm);
  nuis.push_back(ns.astroDeltaGamma);
  nuis.push_back(ns.nuxs);
  nuis.push_back(ns.nubarxs);
  nuis.push_back(ns.kaonLosses);
  nuis.push_back(ns.astroNormSec);
  nuis.push_back(ns.astroDeltaGammaSec);

  assert(nuis.size() == 27);

  return nuis;
}

// Given a human readable flag set, make a bool vector
std::vector<bool> GolemFit::ConvertFitParametersFlag(FitParametersFlag ns) const {
  std::vector<bool> nuis;

  nuis.push_back(ns.convNorm);
  nuis.push_back(ns.promptNorm);
  nuis.push_back(ns.CRDeltaGamma);
  nuis.push_back(ns.piKRatio);
  nuis.push_back(ns.NeutrinoAntineutrinoRatio);
  nuis.push_back(ns.zenithCorrection);

  nuis.push_back(ns.barrWP);
  nuis.push_back(ns.barrWM);
  nuis.push_back(ns.barrYP);
  nuis.push_back(ns.barrYM);
  nuis.push_back(ns.barrZP);
  nuis.push_back(ns.barrZM);
  nuis.push_back(ns.barrHP);
  nuis.push_back(ns.barrHM);

  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.hqdomEfficiency);
  nuis.push_back(ns.holeiceForward);
  nuis.push_back(ns.icegrad0);
  nuis.push_back(ns.icegrad1);
  nuis.push_back(ns.icegrad2);
  nuis.push_back(ns.astroNorm);
  nuis.push_back(ns.astroDeltaGamma);
  nuis.push_back(ns.nuxs);
  nuis.push_back(ns.nubarxs);
  nuis.push_back(ns.kaonLosses);
  nuis.push_back(ns.astroNormSec);
  nuis.push_back(ns.astroDeltaGammaSec);

  assert(nuis.size() == 27);

  return nuis;
}

// And go back to human readable
FitParameters GolemFit::ConvertVecToFitParameters(std::vector<double> vecns) const {

  FitParameters ns;
  ns.convNorm=vecns[0];
  ns.promptNorm=vecns[1];
  ns.CRDeltaGamma=vecns[2];
  ns.piKRatio=vecns[3];
  ns.NeutrinoAntineutrinoRatio=vecns[4];
  ns.zenithCorrection=vecns[5];

  ns.barrWP=vecns[6];
  ns.barrWM=vecns[7];
  ns.barrYP=vecns[8];
  ns.barrYM=vecns[9];
  ns.barrZP=vecns[10];
  ns.barrZM=vecns[11];
  ns.barrHP=vecns[12];
  ns.barrHM=vecns[13];

  ns.domEfficiency=vecns[14];
  ns.hqdomEfficiency=vecns[15];
  ns.holeiceForward=vecns[16];
  ns.icegrad0=vecns[17];
  ns.icegrad1=vecns[18];
  ns.icegrad2=vecns[19];
  ns.astroNorm=vecns[20];
  ns.astroDeltaGamma=vecns[21];
  ns.nuxs=vecns[22];
  ns.nubarxs=vecns[23];
  ns.kaonLosses=vecns[24];
  ns.astroNormSec=vecns[25];
  ns.astroDeltaGammaSec=vecns[26];

  return ns;
}

// Report back the status of object construction
void GolemFit::ReportStatus() const {
  std::cout<< "Data loaded:                   " << CheckDataLoaded() <<std::endl;
  std::cout<< "Sim loaded:                    " << CheckSimulationLoaded()  <<std::endl;
  if(steeringParams_.simType == MCType::LeptonInjector){
    std::cout<< "Dom eff spline constructed:    " << CheckDOMEfficiencySplinesConstructed()<<std::endl;
    std::cout<< "XS weighter constructed:       " << CheckCrossSectionWeighterConstructed() <<std::endl;
    std::cout<< "Flux weighter constructed:     " << CheckFluxWeighterConstructed()<<std::endl;
    std::cout<< "Lepton weighter constructed:   " << CheckLeptonWeighterConstructed()<<std::endl;
    std::cout<< "Oversize weighter constructed: " << CheckOversizeWeighterConstructed()<<std::endl;
  }
  std::cout<< "Data histogram constructed:    " << CheckDataHistogramConstructed()<<std::endl;
  std::cout<< "Sim histogram constructed:     " << CheckSimulationHistogramConstructed()<<std::endl;
  std::cout<< "LLH problem constructed:       " << CheckLikelihoodProblemConstruction()<<std::endl;
}

std::string GolemFit::CheckedFilePath(std::string FilePath) const {
  if(!steeringParams_.quiet) std::cout<<"Reading a file from path "<<FilePath<<std::endl;
  try{
    std::ifstream thefile(FilePath);
    if(thefile.good())
      return FilePath;
    else
      throw std::runtime_error("File " + FilePath + " does not exist!");
  }
  catch(std::exception &re)
    {
      throw std::runtime_error("File " + FilePath + " does not exist!");
    }
}

/*************************************************************************************************************
 * Functions to spit out event distributions and swallow
 * **********************************************************************************************************/

double GolemFit::Swallow(marray<double,2> Data) {
  double TotalWeight=0;
  sample_.clear();
  for(size_t i=0; i!=Data.extent(0); ++i)
    {
      Event e;
      e.energy       = Data[i][0];
      e.zenith       = Data[i][1];
      e.azimuth      = Data[i][2];
      e.topology     = Data[i][3];
      e.year         = Data[i][4];
      //e.length       = Data[i][5];
      e.cachedWeight = Data[i][6];
      TotalWeight+=Data[i][6];
      sample_.push_back(e);
    }
  // remaking data histogram
  if(!steeringParams_.quiet) std::cout<<"Remaking data hist" <<std::endl;
  ConstructDataHistogram();
  // TO DO Improve this
  if(!steeringParams_.quiet) std::cout<<"Reconstrucing likelihood problem" <<std::endl;
  ConstructLikelihoodProblem(priors_, fitSeed_, fixedParams_);
  return TotalWeight;
}

marray<double,2> GolemFit::SpitData() const {
  marray<double,2> ReturnVec { sample_.size(), 7} ;
  for(size_t i=0; i!=sample_.size(); ++i) {
      ReturnVec[i][0]=sample_[i].energy;
      ReturnVec[i][1]=sample_[i].zenith;
      ReturnVec[i][2]=sample_[i].azimuth;
      ReturnVec[i][3]=sample_[i].topology;
      ReturnVec[i][4]=sample_[i].year;
      //ReturnVec[i][5]=sample_[i].length;
      ReturnVec[i][5]=1.0;
      ReturnVec[i][6]=sample_[i].cachedWeight;
  }
  return ReturnVec;
}

double GolemFit::SetupDataChallenge(int seed, FitParameters fit_params_dc) {
  // Save this thing for later, since we'll be adjusting it
  bool my_quiet=steeringParams_.quiet;

  // Set the parameter point to the data challenge point and make realization
  steeringParams_.quiet=true;
  marray<double,2> TheRealization=SpitRealization(ConvertFitParameters(fit_params_dc),seed);

  // Set the parameter point back to the original one
  steeringParams_.quiet=my_quiet;

  // Set the realization as data
  return Swallow(TheRealization);
}

double GolemFit::SetupDataChallenge(int seed, FitParameters fit_params_dc, NewPhysicsParams npp_dc) {
  // Save these two things for later, since we'll be adjusting them
  bool my_quiet=steeringParams_.quiet;
  NewPhysicsParams my_npp=new_physics_params_;

  // Set the parameter point to the data challenge point and make realization
  steeringParams_.quiet=true;
  SetNewPhysicsParams(npp_dc);
  marray<double,2> TheRealization=SpitRealization(ConvertFitParameters(fit_params_dc),seed);

  // Set the parameter point back to the original one
  SetNewPhysicsParams(my_npp);
  steeringParams_.quiet=my_quiet;

  // Set the realization as data
  return Swallow(TheRealization);
}

marray<double,2> GolemFit::SpitRealization( FitParameters nuisance, int seed) const {
  return SpitRealization(ConvertFitParameters(nuisance), seed);
}

marray<double,2> GolemFit::SpitRealization( std::vector<double> nuisance, int seed) const {
  std::mt19937 rng;
  rng.seed(seed);

  auto weighter=DFWM(nuisance);
  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    if(std::isnan(w) || std::isinf(w) || w<0){
      std::cout << ConvertVecToFitParameters(nuisance) << std::endl;
      std::cout << "Bad weight! " << w << std::endl;
      std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
      std::cout << e.energy << ' ' << e.year << ' ' << w << std::endl;
    }
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization=phys_tools::likelihood::generateSample(weights,mainSimulation_,expected,rng);

  marray<double,2> ReturnVec { realization.size(), 7} ;

  for(size_t i=0; i!=realization.size(); ++i)
    {
      ReturnVec[i][0]=realization[i].energy;
      ReturnVec[i][1]=realization[i].zenith;
      ReturnVec[i][2]=realization[i].azimuth;
      ReturnVec[i][3]=realization[i].topology;
      ReturnVec[i][4]=realization[i].year;
      //ReturnVec[i][5]=realization[i].length;
      ReturnVec[i][5]=1.;
      ReturnVec[i][6]=1.;
    }
  return ReturnVec;
}

marray<double,2> GolemFit::SpitExpectation(FitParameters fit_params) const {
  return SpitExpectation(ConvertFitParameters(fit_params));
}

marray<double,2> GolemFit::SpitExpectation(std::vector<double> fit_params) const {
  const std::deque<Event> & sim = (fastmode_constructed_)?metaEvents_:mainSimulation_;
  marray<double,2> ReturnVec {sim.size(), 7} ;
  auto weighter = DFWM(fit_params);
  for(size_t i=0; i!=sim.size(); ++i) {
      ReturnVec[i][0]=sim[i].energy;
      ReturnVec[i][1]=sim[i].zenith;
      ReturnVec[i][2]=sim[i].azimuth;
      ReturnVec[i][3]=sim[i].topology;
      ReturnVec[i][4]=sim[i].year;
      //ReturnVec[i][5]=sim[i].length;
      ReturnVec[i][5]=1.;
      ReturnVec[i][6]=weighter(sim[i]);
  }
  return ReturnVec;
}

std::deque<Event> GolemFit::SpitMonteCarlo() const {
  return (fastmode_constructed_)?metaEvents_:mainSimulation_;
}

void GolemFit::SpitWeightedMonteCarloToHDF5(FitParameters fit_parameters, std::string str) const {
  hid_t file_id,root_id;
  hid_t dset_id;

  if(not steeringParams_.quiet)
    std::cout << "Beginning spitting in HDF5 " << str << std::endl;

  // create file
  file_id = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file_id < 0)
    throw std::runtime_error("GolemFit::Error::Cannot create file at " + str + ".");
  root_id = H5Gopen(file_id, "/",H5P_DEFAULT);

  auto serialize_to_hdf5 = [&](std::function<double(const Event&)> get_variable, std::string dataset_name){
    std::vector<double> tmp_vector(mainSimulation_.size());
    for(size_t i=0; i < mainSimulation_.size(); i++)
      tmp_vector[i] = get_variable(mainSimulation_[i]);
    hsize_t dim[1] {static_cast<hsize_t>(mainSimulation_.size())};
    dset_id = H5LTmake_dataset(root_id,dataset_name.c_str(),1,dim,H5T_NATIVE_DOUBLE,static_cast<const void*>(tmp_vector.data()));
  };

  auto serialize_vector_to_hdf5 = [&](const std::vector<double> & array, std::string dataset_name){
    hsize_t dim[1] {static_cast<hsize_t>(array.size())};
    dset_id = H5LTmake_dataset(root_id,dataset_name.c_str(),1,dim,H5T_NATIVE_DOUBLE,static_cast<const void*>(array.data()));
  };

  // event unique identifier
  serialize_vector_to_hdf5(run_number_,"run");
  serialize_vector_to_hdf5(event_number_,"event");

  serialize_vector_to_hdf5(subevent_number_,"subevent");
  // true event properties
  serialize_to_hdf5([&](const Event & e)->double {return (int)e.primaryType;},"primaryType");
  serialize_to_hdf5([&](const Event & e)->double {return e.primaryEnergy;},"primaryEnergy");
  serialize_to_hdf5([&](const Event & e)->double {return e.primaryAzimuth;},"primaryAzimuth");
  serialize_to_hdf5([&](const Event & e)->double {return e.primaryZenith;},"primaryZenith");

  // observable event properties
  serialize_to_hdf5([&](const Event & e)->double {return e.energy;}, "energy");
  serialize_to_hdf5([&](const Event & e)->double {return e.zenith;}, "zenith");
  serialize_to_hdf5([&](const Event & e)->double {return e.azimuth;}, "azimuth");

  // all the weights
  serialize_to_hdf5([&](const Event & e)->double {return e.oneWeight;},"oneWeight");
  serialize_to_hdf5([&](const Event & e)->double {return e.number_of_generated_mc_events;},"number_of_generated_mc_events");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedLivetime;},"cachedLivetime");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedWeight;},"cachedWeight");

  serialize_to_hdf5([&](const Event & e)->double {return e.cachedConvPionWeight;},"cachedConvPionWeight");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedConvKaonWeight;},"cachedConvKaonWeight");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedConvWeight;},"cachedConvWeight");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedPromptWeight;},"cachedPromptWeight");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedAstroMuWeight;},"cachedAstroMuWeight");

  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModWP;},"cachedBarrModWP");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModWM;},"cachedBarrModWM");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModYP;},"cachedBarrModYP");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModYM;},"cachedBarrModYM");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModZP;},"cachedBarrModZP");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedBarrModZM;},"cachedBarrModZM");

  serialize_to_hdf5([&](const Event & e)->double {return e.cachedHoleIceConv;},"cachedHoleIceConv");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedHoleIcePrompt;},"cachedHoleIcePrompt");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedHoleIceAstro;},"cachedHoleIceAstro");

  serialize_to_hdf5([&](const Event & e)->double {return e.cachedDOMEffConv;},"cachedDOMEffConv");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedDOMEffPrompt;},"cachedDOMEffPrompt");
  serialize_to_hdf5([&](const Event & e)->double {return e.cachedDOMEffAstro;},"cachedDOMEffAstro");

  auto weighter = DFWM(ConvertFitParameters(fit_parameters));
  serialize_to_hdf5([&](const Event & e)->double {return weighter(e);},"Weight");

  // close root group
  H5Gclose(root_id);
  // close HDF5 file
  H5Fclose(file_id);
  // close all HDF5 variables/memory
  H5close();

  if(not steeringParams_.quiet)
    std::cout << "End spitting in HDF5 " << str << std::endl;
}

bool GolemFit::SetupAsimov(FitParameters fit_parameters) {
  if(!steeringParams_.quiet){
    std::cout<<"Setting up Asimov set with parameters" <<std::endl;
    std::cout << fit_parameters << std::endl;
  }
  SetupAsimov(ConvertFitParameters(fit_parameters));
  return true;
}

void GolemFit::SetupAsimov(std::vector<double> fit_parameters) {
  Swallow(SpitExpectation(fit_parameters));
  asimovSeed_ = ConvertVecToFitParameters(fit_parameters);
//  asimov_setup_ = true;
}

bool GolemFit::CheckAsimovSanity() const{
  return !(asimov_setup_ && steeringParams_.allow_mismatched_parameters && !check_fit_parameters_equality(fitSeed_, asimovSeed_, fixedParams_));
}

/*************************************************************************************************************
 * Functions to get bin edges
 * **********************************************************************************************************/

// Given a histogram reference h, get bin edges in dimension dim
 template<unsigned int hist_index>
 std::vector<double> GolemFit::PullBinEdges(int dim, const HistType& hh) const{
   auto h = std::get<hist_index>(hh);
   std::vector<double> edges_i(h.getBinCount(dim));
   for(unsigned int j=0; j<h.getBinCount(dim); j++)
     edges_i[j]=h.getBinEdge(dim,j);
   edges_i.push_back(h.getBinEdge(dim,h.getBinCount(dim)-1)+h.getBinWidth(dim,h.getBinCount(dim)-1));
   return edges_i;
 }

 std::vector<double> GolemFit::GetEnergyBinsData() const{
   return PullBinEdges<0>(0,dataHist_);
 }
 std::vector<double> GolemFit::GetEnergyBinsMC() const{
   return PullBinEdges<0>(0,simHist_);
 }

 std::vector<double> GolemFit::GetZenithBinsData() const{
   return PullBinEdges<0>(1,dataHist_);
 }
 std::vector<double> GolemFit::GetZenithBinsMC() const{
   return PullBinEdges<0>(1,simHist_);
 }

 std::vector<double> GolemFit::GetAzimuthBinsData() const{
   return PullBinEdges<0>(2,dataHist_);
 }
 std::vector<double> GolemFit::GetAzimuthBinsMC() const{
   return PullBinEdges<0>(2,simHist_);
 }

 std::vector<double> GolemFit::GetTopologyBinsData() const{
   return PullBinEdges<0>(3,dataHist_);
 }
 std::vector<double> GolemFit::GetTopologyBinsMC() const{
   return PullBinEdges<0>(3,simHist_);
 }

 std::vector<double> GolemFit::GetYearBinsData() const{
     return PullBinEdges<0>(4,dataHist_);
 }
 std::vector<double> GolemFit::GetYearBinsMC() const{
     return PullBinEdges<0>(4,simHist_);
 }

 std::vector<std::vector<std::vector<double>>> GolemFit::GetBinsData() const{
   std::vector<std::vector<std::vector<double>>> all_bins(2);

   std::vector<std::vector<double>>& ct_bins = all_bins[0];
   std::vector<std::vector<double>>& aux_bins = all_bins[1];

   ct_bins.reserve(5);
   aux_bins.reserve(5);

   ct_bins[0] = PullBinEdges<0>(0, dataHist_);
   ct_bins[1] = PullBinEdges<0>(1, dataHist_);
   ct_bins[2] = PullBinEdges<0>(2, dataHist_);
   ct_bins[3] = PullBinEdges<0>(3, dataHist_);
   ct_bins[4] = PullBinEdges<0>(4, dataHist_);

   aux_bins[0] = PullBinEdges<0>(0, auxdataHist_);
   aux_bins[1] = PullBinEdges<0>(1, auxdataHist_);
   aux_bins[2] = PullBinEdges<0>(2, auxdataHist_);
   aux_bins[3] = PullBinEdges<0>(3, auxdataHist_);
   aux_bins[4] = PullBinEdges<0>(4, auxdataHist_);

   return all_bins;
 }

 std::vector<std::vector<std::vector<double>>> GolemFit::GetBinsMC() const{
   std::vector<std::vector<std::vector<double>>> all_bins(2);

   std::vector<std::vector<double>>& ct_bins = all_bins[0];
   std::vector<std::vector<double>>& aux_bins = all_bins[1];

   ct_bins.reserve(5);
   aux_bins.reserve(5);

   ct_bins[0] = PullBinEdges<0>(0, simHist_);
   ct_bins[1] = PullBinEdges<0>(1, simHist_);
   ct_bins[2] = PullBinEdges<0>(2, simHist_);
   ct_bins[3] = PullBinEdges<0>(3, simHist_);
   ct_bins[4] = PullBinEdges<0>(4, simHist_);

   aux_bins[0] = PullBinEdges<0>(0, auxsimHist_);
   aux_bins[1] = PullBinEdges<0>(1, auxsimHist_);
   aux_bins[2] = PullBinEdges<0>(2, auxsimHist_);
   aux_bins[3] = PullBinEdges<0>(3, auxsimHist_);
   aux_bins[4] = PullBinEdges<0>(4, auxsimHist_);

   return all_bins;
 }

/*************************************************************************************************************
 * Functions to construct nusquids state on the fly
 * **********************************************************************************************************/

void GolemFit::ConstructNuSQuIDSObjects(){
  bool interaction = (steeringParams_.simType == MCType::LeptonInjector);
  switch(new_physics_params_.type){
    case NewPhysicsType::None: ConstructNuSQuIDSObjectsForStandardModel(interaction); return;
    case NewPhysicsType::NonStandardInteraction: ConstructNuSQuIDSObjectsForNSI(interaction); return;
    case NewPhysicsType::LorentzViolation: ConstructNuSQuIDSObjectsForLV(interaction); return;
    case NewPhysicsType::NuSQuIDS: ConstructNuSQuIDSObjectsForCrossSectionModification(true); return;
    default: if(!steeringParams_.quiet) std::cout << "No nuSQuIDS calculation performed/needed for set NewPhysicsType" << std::endl; return;
  };
}

std::map<FluxComponent,nusquids::marray<double,4>> GolemFit::ConstructNuSQuIDSFluxComponentArrays(nusquids::marray<double,1> costh_nodes, nusquids::marray<double,1> energy_nodes) const {
  marray<double,4> inistate_kaon   {costh_nodes.size(),energy_nodes.size(),2,numneu};
  marray<double,4> inistate_pion   {costh_nodes.size(),energy_nodes.size(),2,numneu};
  marray<double,4> inistate_prompt {costh_nodes.size(),energy_nodes.size(),2,numneu};
  marray<double,4> inistate_astro  {costh_nodes.size(),energy_nodes.size(),2,numneu};

  LW::Event scratch_lw_e;
  for(unsigned  int ci = 0 ; ci < costh_nodes.size(); ci++){
    for(unsigned  int ei = 0 ; ei < energy_nodes.size(); ei++){
      double enu = energy_nodes[ei]/units.GeV;
      double cth = costh_nodes[ci];

      scratch_lw_e.energy=enu;
      scratch_lw_e.zenith=acos(cth);

      scratch_lw_e.primary_type=LW::ParticleType::NuE;
      inistate_kaon[ci][ei][0][0] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][0][0] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][0][0] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][0][0] = (*fluxAstro_)(scratch_lw_e)*sclup_;
      scratch_lw_e.primary_type=LW::ParticleType::NuMu;
      inistate_kaon[ci][ei][0][1] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][0][1] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][0][1] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][0][1] = (*fluxAstro_)(scratch_lw_e)*sclup_;
      scratch_lw_e.primary_type=LW::ParticleType::NuTau;
      inistate_kaon[ci][ei][0][2] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][0][2] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][0][2] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][0][2] = (*fluxAstro_)(scratch_lw_e)*sclup_;

      scratch_lw_e.primary_type=LW::ParticleType::NuEBar;
      inistate_kaon[ci][ei][1][0] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][1][0] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][1][0] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][1][0] = (*fluxAstro_)(scratch_lw_e)*sclup_;
      scratch_lw_e.primary_type=LW::ParticleType::NuMuBar;
      inistate_kaon[ci][ei][1][1] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][1][1] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][1][1] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][1][1] = (*fluxAstro_)(scratch_lw_e)*sclup_;
      scratch_lw_e.primary_type=LW::ParticleType::NuTauBar;
      inistate_kaon[ci][ei][1][2] = (*fluxKaon_)(scratch_lw_e)*sclup_;
      inistate_pion[ci][ei][1][2] = (*fluxPion_)(scratch_lw_e)*sclup_;
      inistate_prompt[ci][ei][1][2] = (*fluxPrompt_)(scratch_lw_e)*sclup_;
      inistate_astro[ci][ei][1][2] = (*fluxAstro_)(scratch_lw_e)*sclup_;
    }
  }

  return {{FluxComponent::atmPion,inistate_pion},{FluxComponent::atmKaon,inistate_kaon},{FluxComponent::atmPrompt,inistate_prompt},{FluxComponent::diffuseAstro,inistate_astro}};
}

void GolemFit::ConstructNuSQuIDSObjectsForNSI(bool interactions){
  // to do shivesh check that things make sense
  lw_nusquids_fluxes_nsi_.clear();
  nusquids::marray<double,1> cos_range = linspace(-1.,1.,50);
  nusquids::marray<double,1> e_range = logspace(1.e4*units.GeV,1.e8*units.GeV,101);
  auto initial_fluxes = ConstructNuSQuIDSFluxComponentArrays(cos_range,e_range);
  auto xs = std::make_shared<tools::ScaledNeutrinoCrossSections>(XSECTION_LOCATION "csms.h5");
  for(auto component : {FluxComponent::atmPion, FluxComponent::atmKaon, FluxComponent::atmPrompt, FluxComponent::diffuseAstro}){
    nusquids::nuSQUIDSAtm<nuSQUIDSNSI> nus_atm_nsi(cos_range,new_physics_params_.epsilon_mutau,e_range,numneu,both,interactions);
    ReconfigureNuSQuIDSObjectSettings(nus_atm_nsi);
    nus_atm_nsi.Set_initial_state(initial_fluxes[component],nusquids::flavor);
    nus_atm_nsi.EvolveState();
    lw_nusquids_fluxes_nsi_.insert({component,std::make_shared<LW::nuSQUIDSAtmFlux<nuSQUIDSNSI>>(std::move(nus_atm_nsi))});
  }
  nusquids_nsi_fluxes_calculated_ = true;
}

void GolemFit::ConstructNuSQuIDSObjectsForStandardModel(bool interactions){
  lw_nusquids_fluxes_.clear();
  nusquids::marray<double,1> cos_range = linspace(-1.,1.,50);
  nusquids::marray<double,1> e_range = logspace(1.e4*units.GeV,1.e8*units.GeV,101);
  auto initial_fluxes = ConstructNuSQuIDSFluxComponentArrays(cos_range,e_range);
  auto xs = std::make_shared<tools::ScaledNeutrinoCrossSections>(XSECTION_LOCATION "csms.h5");
  for(auto component : {FluxComponent::atmPion, FluxComponent::atmKaon, FluxComponent::atmPrompt, FluxComponent::diffuseAstro}){
    nusquids::nuSQUIDSAtm<> nus_atm(cos_range,e_range,numneu,both,interactions,xs);
    ReconfigureNuSQuIDSObjectSettings(nus_atm, false);
    nus_atm.Set_initial_state(initial_fluxes[component],nusquids::flavor);
    nus_atm.EvolveState();
    lw_nusquids_fluxes_.insert({component,std::make_shared<LW::nuSQUIDSAtmFlux<>>(std::move(nus_atm))});
  }
  nusquids_nsi_fluxes_calculated_ = true;
}

void GolemFit::ConstructNuSQuIDSObjectsForLV(bool interactions){
  // to do shivesh check that things make sense
  if (!steeringParams_.quiet) std::cout << "Entering ConstructNuSQuIDSObjectsForLV." << std::endl;
  lw_nusquids_fluxes_lv_.clear();
  nusquids::marray<double,1> cos_range = linspace(-1.,1.,50);
  nusquids::marray<double,1> e_range = logspace(1.e4*units.GeV,1.e8*units.GeV,101);
  auto initial_fluxes = ConstructNuSQuIDSFluxComponentArrays(cos_range,e_range);
  auto xs = std::make_shared<tools::ScaledNeutrinoCrossSections>(XSECTION_LOCATION "csms.h5");

  // calculate lorentz violating operator
  squids::Const squids_parameter_container;
  // make the inner diagonal piece of the operator
  squids::SU_vector LVOp = (1./(new_physics_params_.lambda_1*pow(units.GeV,new_physics_params_.n_lv-1.)))*squids::SU_vector::Projector(numneu,1); // diag(0,1/lambda1,0)
  LVOp += (1./new_physics_params_.lambda_2)*squids::SU_vector::Projector(numneu,2);// diag(0,1/lambda1,1/lambda2)
  // set the rotations to be applied to the operator
  squids_parameter_container.SetMixingAngle(0,1,new_physics_params_.th12);
  squids_parameter_container.SetMixingAngle(0,2,new_physics_params_.th13);
  squids_parameter_container.SetMixingAngle(1,2,new_physics_params_.th23);
  // rotate it by the unitary transformation given by those angles
  LVOp.RotateToB1(squids_parameter_container);
  if (!steeringParams_.quiet) std::cout << "Rotated Unitary Transform" << std::endl;

  for(auto component : {FluxComponent::atmPion, FluxComponent::atmKaon, FluxComponent::atmPrompt, FluxComponent::diffuseAstro}){
    if (!steeringParams_.quiet) std::cout << "In for loop" << std::endl;
    nusquids::nuSQUIDSAtm<nuSQUIDSLV> nus_atm_lv(cos_range,e_range,numneu,both,interactions);
    ReconfigureNuSQuIDSObjectSettings(nus_atm_lv);
    nus_atm_lv.Set_initial_state(initial_fluxes[component],nusquids::flavor);
    // the LV class is not as nice as the NSI one since you have to set
    // a bunch of things manually; maybeb want to improve the LV.h constructor like the NSI one
    // this loop access the internal nusquids elements and does the magic. CA
    std::vector<nuSQUIDSLV>& nus_atm_lv_inner_array = nus_atm_lv.GetnuSQuIDS();
    if (!steeringParams_.quiet) std::cout << "Entering inner for loop" << std::endl;
    for(nuSQUIDSLV& nus_atm_lv_inner : nus_atm_lv_inner_array){
      nus_atm_lv_inner.Set_LV_EnergyPower(new_physics_params_.n_lv);
      nus_atm_lv_inner.Set_LV_Operator(LVOp);
    }
    if (!steeringParams_.quiet) std::cout << "Leaving inner for loop" << std::endl;
    nus_atm_lv.EvolveState();
    lw_nusquids_fluxes_lv_.insert({component,std::make_shared<LW::nuSQUIDSAtmFlux<nuSQUIDSLV>>(std::move(nus_atm_lv))});
    if (!steeringParams_.quiet) std::cout << "Doing next iteration" << std::endl;
  }
  if (!steeringParams_.quiet) std::cout << "Leaving for loop" << std::endl;
  nusquids_lv_fluxes_calculated_ = true;
  if (!steeringParams_.quiet) std::cout << "Leaving ConstructNuSQuIDSObjectsForLV." << std::endl;
}

void GolemFit::ConstructNuSQuIDSObjectsForCrossSectionModification(bool interactions){
  lw_nusquids_fluxes_xs_.clear();
  nusquids::marray<double,1> cos_range = linspace(-1.,1.,50);
  nusquids::marray<double,1> e_range = logspace(1.e4*units.GeV,1.e8*units.GeV,101);
  auto initial_fluxes = ConstructNuSQuIDSFluxComponentArrays(cos_range,e_range);
  // reinitialize nusquidsatm objects because std::move will
  // invalidate them for later use
  // setup modified cross section stuff
  std::vector<double> edges = {60*units.TeV, 100*units.TeV, 200*units.TeV, 500*units.TeV, 10000*units.TeV};
  std::vector<double> scales = {new_physics_params_.xs0, new_physics_params_.xs1, new_physics_params_.xs2, new_physics_params_.xs3};  
  xsv = std::make_shared<tools::ScaledNeutrinoCrossSections>(XSECTION_LOCATION "csms.h5", edges, scales);
  for(auto component : {FluxComponent::atmPion, FluxComponent::atmKaon, FluxComponent::atmPrompt, FluxComponent::diffuseAstro}){
    nusquids::nuSQUIDSAtm<> nus_atm(cos_range,e_range,numneu,both,interactions,xsv);
    ReconfigureNuSQuIDSObjectSettings(nus_atm, false);
    nus_atm.Set_initial_state(initial_fluxes[component],nusquids::flavor);
    nus_atm.EvolveState();
    lw_nusquids_fluxes_xs_.insert({component,std::make_shared<LW::nuSQUIDSAtmFlux<>>(std::move(nus_atm))});
  }
}

std::tuple<double,double> GolemFit::EquatorialToGalacticCoordinates(double ra, double dec) const {
  if(std::isnan(ra) or std::isnan(dec))
    throw std::runtime_error("Directions should be numbers.");
  // returns galactic longitude (l) and latitude (b) in epoch J2000;
  double radians = M_PI/180.;
  double alpha = 192.859508*radians;
  double delta = 27.128336*radians;
  double la = 33.0*radians;

  double b = asin(sin(dec)*sin(delta) + cos(dec)*cos(delta)*cos(ra-alpha));
  double l = atan2(sin(dec)*cos(delta) - cos(dec) * sin(delta) * cos(ra-alpha),
                   cos(dec)*sin(ra-alpha)) + la;

  double rl = std::remainder(l,2.*M_PI);
  return std::make_tuple(rl, b);
}

} // close namespace golemfit
