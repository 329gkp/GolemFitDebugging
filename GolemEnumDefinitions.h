#ifndef GOLEMENUMDEFINITIONS_H_
#define GOLEMENUMDEFINITIONS_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <LeptonWeighter/ParticleType.h>

namespace golemfit{

enum class CompactMode {dump, jason};

enum class MCType {LeptonInjector, NuGen, MuonGun};

enum class sampleTag{HESE, MESE, Sterile, MagicTau};
std::string GetSampleName(sampleTag sample);
std::ostream& operator<<(std::ostream& os, const sampleTag& sample);

enum class AstrophysicalNeutrinoModel {SteckerAGN, Fang_UHECR, LLGAGN_modelB1, LLGAGN_modelB4, LLGAGN_twocompModel, Maria_BLLacs, Murase_chockedJets, SBG_minBmodel, Tavecchi_lowPower, TDE_WinterBiehl, Globus_GZK_GRBEvol_MixComp_beta2, Globus_GZK_GRBEvol_MixComp_beta25 , Globus_GZK_GRBEvol_Proton, Globus_GZK_NoEvol_MixComp_beta25, Globus_GZK_NoEvol_MixComp_beta2, Globus_GZK_NoEvol_Proton, Ahlers_GZK};
std::string GetAstrophysicalNeutrinoModelName(AstrophysicalNeutrinoModel model);
std::ostream& operator<<(std::ostream& os, const AstrophysicalNeutrinoModel& model);

enum class EarthPropagationInterfaces {nuSQuIDS,nuFATE};

enum class BarrParameter {GP,GM,IP,IM,HP,HM,WP,WM,YP,YM,ZP,ZM};
std::string GetBarrParameterName(const BarrParameter& barr_parameter);
std::ostream& operator<<(std::ostream& os, const BarrParameter& barr_parameter);

enum class FluxComponent {atmMuon, atmPion, atmKaon, atmConv, atmPrompt, diffuseAstro, galacticAstro, diffuseAstro_e, diffuseAstro_mu, diffuseAstro_tau};
std::string GetFluxComponentName(const FluxComponent& component);
std::ostream& operator<<(std::ostream& os, const FluxComponent& component);

enum class Topology {shower = 0, track = 1, doublebang = 2};
std::string GetTopologyName(const Topology& topology);
std::ostream& operator<<(std::ostream& os, const Topology& model);

enum class DiffuseFitType {SinglePowerLaw, DoublePowerLaw, LogParaboloid, AdHocModel, Galactic, SinglePowerDoubleNorm, DoublePowerSingleNorm};

enum class NewPhysicsType {None, DarkMatterAnnihilation, DarkMatterScattering,
                           DarkMatterDecay, LorentzViolation, NonStandardInteraction,
                           SterileNeutrino, NewCrossSection, NuSQuIDS};

std::string GetNewPhysicsTypeName(const NewPhysicsType&);

// this function should totally get moved to LeptonWeighter CA. Note this aschneider on the review thingy.
std::string GetParticleName(const LW::ParticleType& particle);

/* don't forget me carlos; think.
 * sorry. I forgot you function. not sure what this was fore. ill think.
struct NuSQuIDSFluxComponentArrays {
  std::map<FluxComponent,nusquids::marray<double,4>> nusquids_initial_flux_array;
};
*/

}

#endif // close namespace
