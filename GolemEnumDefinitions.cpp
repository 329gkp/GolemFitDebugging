#include "GolemEnumDefinitions.h"

namespace golemfit {

std::string GetAstrophysicalNeutrinoModelName(AstrophysicalNeutrinoModel model){
  switch(model) {
    case AstrophysicalNeutrinoModel::SteckerAGN : return "stecker_agn";
    case AstrophysicalNeutrinoModel::Fang_UHECR : return "fang_uhecr";
    case AstrophysicalNeutrinoModel::LLGAGN_modelB1 : return "LLGAGN_modelB1";
    case AstrophysicalNeutrinoModel::LLGAGN_modelB4 : return "LLGAGN_modelB4";
    case AstrophysicalNeutrinoModel::LLGAGN_twocompModel : return "LLGAGN_twocompModel";
    case AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_MixComp_beta2: return "Globus_GZK_NoEvol_Mixcomp_beta2";
    case AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_MixComp_beta25: return "Globus_GZK_NoEvol_Mixcomp_beta25";
    case AstrophysicalNeutrinoModel::Globus_GZK_NoEvol_Proton: return "Globus_GZK_NoEvol_proton";
    case AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_MixComp_beta2: return "Globus_GZK_GRBevol_Mixcomp_beta2";
    case AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_MixComp_beta25: return "Globus_GZK_GRBevol_Mixcomp_beta25";
    case AstrophysicalNeutrinoModel::Globus_GZK_GRBEvol_Proton: return "Globus_GZK_GRBEvol_proton";
    case AstrophysicalNeutrinoModel::Maria_BLLacs: return "maria_BLLacs";
    case AstrophysicalNeutrinoModel::Murase_chockedJets: return "Murase+_chokedJets";
    case AstrophysicalNeutrinoModel::SBG_minBmodel: return "SBG_minBmodel";
    case AstrophysicalNeutrinoModel::Tavecchi_lowPower: return "Tavecchi_lowPower";
    case AstrophysicalNeutrinoModel::TDE_WinterBiehl: return "TDE_WinterBiehl";
    case AstrophysicalNeutrinoModel::Ahlers_GZK: return "Ahlers_GZK";
    default: return std::string();
  }
}

std::ostream& operator<<(std::ostream& os, const AstrophysicalNeutrinoModel& model){
  os << GetAstrophysicalNeutrinoModelName(model);
  return os;
}

std::string GetSampleName(sampleTag sample){
  switch(sample) {
    case sampleTag::HESE: return "HESE";
    case sampleTag::MESE: return "MESE";
    case sampleTag::Sterile: return "STERILE";
    case sampleTag::MagicTau: return "MAGICTAU";
    default: return std::string();
  }
}

std::ostream& operator<<(std::ostream& os, const sampleTag& sample){
  os << GetSampleName(sample);
  return os;
}

std::string GetBarrParameterName(const BarrParameter& barr_parameter){
  switch(barr_parameter) {
    // mainly pion paramters
    case BarrParameter::GP: return "GP";
    case BarrParameter::GM: return "GM";
    case BarrParameter::IP: return "IP";
    case BarrParameter::IM: return "IM";
    case BarrParameter::HP: return "HP";
    case BarrParameter::HM: return "HM";
    // mainly kaon paramters
    case BarrParameter::WP: return "WP";
    case BarrParameter::WM: return "WM";
    case BarrParameter::YP: return "YP";
    case BarrParameter::YM: return "YM";
    case BarrParameter::ZP: return "ZP";
    case BarrParameter::ZM: return "ZM";
    default: throw std::runtime_error("Barr paramter given by enumerator not registered. No label guy found.");
  }
}

std::ostream& operator<<(std::ostream& os, const BarrParameter& barr_parameter){
  os << GetBarrParameterName(barr_parameter);
  return os;
}

std::string GetFluxComponentName(const FluxComponent& component){
  switch(component) {
    case FluxComponent::atmMuon: return "atmMuon";
    case FluxComponent::atmPion: return "atmPion";
    case FluxComponent::atmKaon: return "atmKaon";
    case FluxComponent::atmConv: return "atmConv";
    case FluxComponent::atmPrompt: return "atmPrompt";
    case FluxComponent::diffuseAstro: return "diffuseAstro";
    case FluxComponent::galacticAstro: return "galacticAstro";
    case FluxComponent::diffuseAstro_e: return "diffuseAstro_e";
    case FluxComponent::diffuseAstro_mu: return "diffuseAstro_mu";
    case FluxComponent::diffuseAstro_tau: return "diffuseAstro_tau";
    default: throw std::runtime_error("Flux component given by enumerator not registered. No label guy found.");
  }
}

std::ostream& operator<<(std::ostream& os, const FluxComponent& component){
  os << GetFluxComponentName(component);
  return os;
}

std::string GetTopologyName(const Topology& topology){
  switch(topology) {
    case Topology::shower: return "shower";
    case Topology::track: return "track";
    case Topology::doublebang: return "doublebang";
    default: return std::string();
  }
}

std::ostream& operator<<(std::ostream& os, const Topology& topology){
  os << GetTopologyName(topology);
  return os;
}

std::string GetNewPhysicsTypeName(const NewPhysicsType& np) {
  switch(np) {
    case NewPhysicsType::None: return "None";
    case NewPhysicsType::DarkMatterAnnihilation: return "DarkMatterAnnihilation";
    case NewPhysicsType::DarkMatterScattering: return "DarkMatterScattering";
    case NewPhysicsType::DarkMatterDecay: return "DarkMatterDecay";
    case NewPhysicsType::LorentzViolation: return "LorentzViolation";
    case NewPhysicsType::NonStandardInteraction: return "NonStandardInteraction";
    case NewPhysicsType::SterileNeutrino: return "SterileNeutrino";
    case NewPhysicsType::NewCrossSection: return "NewCrossSection";
    case NewPhysicsType::NuSQuIDS: return "NuSQuIDS";
    default: return std::string();
  }
}

std::string GetParticleName(const LW::ParticleType& particle){
  switch(particle) {
    case LW::ParticleType::NuE: return "nue";
    case LW::ParticleType::NuEBar: return "nuebar";
    case LW::ParticleType::NuMu: return "numu";
    case LW::ParticleType::NuMuBar: return "numubar";
    case LW::ParticleType::NuTau: return "nutau";
    case LW::ParticleType::NuTauBar: return "nutaubar";
    default: return std::string();
  }
}

} // close namespace
