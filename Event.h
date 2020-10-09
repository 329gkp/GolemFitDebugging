#ifndef EVENT_H
#define EVENT_H

#include <iostream>
#include <set>

#include <boost/math/constants/constants.hpp>
#include <LeptonWeighter/ParticleType.h>
#include <PhysTools/tableio.h>
#include "GolemEnumDefinitions.h"
#include "FastMode.h"
#include "json.hpp"

double differenceAngle(double zenith1, double azimuth1, double zenith2, double azimuth2);
herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data);

struct vector{
	double x,y,z;
	vector():x(0),y(0),z(0){}
	vector(double x, double y, double z):x(x),y(y),z(z){}
	vector(double zenith, double azimuth){
		using boost::math::constants::pi;
		const double theta = pi<double>()-zenith;
		const double phi = azimuth-pi<double>();
		const double rho = std::sin(theta);
		x = rho*std::cos(phi);
		y = rho*std::sin(phi);
		z = std::cos(theta);
	}

	double operator*(const vector& rhs) const{
		return x*rhs.x + y*rhs.y + z*rhs.z;
	}
	vector& operator*=(double a){
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	vector operator*(double a) const{
		return vector(*this)*=a;
	}
	vector& operator-=(const vector& rhs){
		x-=rhs.x;
		y-=rhs.y;
		z-=rhs.z;
		return *this;
	}
	vector operator-(const vector& rhs) const{
		return vector(*this)-=rhs;
	}
	double magnitude() const{
		return sqrt(*this * *this);
	}
};

enum class Level {generation,level_1,level_2,neutrino};

class babyEvent {
  public:
  // LW/LI compatibility
  LW::ParticleType primaryType;
  float primaryEnergy;
  float primaryAzimuth;
  float primaryZenith;
  // reconstructed quantities
  float energy;
  float zenith;
  float azimuth;
  // weight
  float oneWeight;
  // topology
  unsigned int topology;
  babyEvent(){}
};

class Event{
public:
  // event identifier
  //unsigned int run;
  //unsigned int event;
  //unsigned int subevent;

  // LW/LI compatibility
  LW::ParticleType final_state_particle_0;
  LW::ParticleType final_state_particle_1;
  LW::ParticleType primaryType;
  float primaryEnergy;
  float primaryAzimuth;
  float primaryZenith;
  float totalColumnDepth;
  float intX,intY;

  // for weighting
  float oneWeight;
  float number_of_generated_mc_events;

  // meta-event quantities
  int num_events;

  // reconstructed quantities
  float energy;
  float zenith;
  float azimuth;

  // which sample
  golemfit::sampleTag sample;

  // topology
  unsigned int topology;

  // cached quantities
  float cachedLivetime;
  float cachedWeight;

  // cached fluxes
  double cachedConvPionWeight;
  double cachedConvKaonWeight;
  double cachedConvWeight;
  double cachedPromptWeight;
  double cachedAstroMuWeight;

  // barr parameterers caches
  double cachedBarrModWP;
  double cachedBarrModWM;
  double cachedBarrModYP;
  double cachedBarrModYM;
  double cachedBarrModZP;
  double cachedBarrModZM;

  // cached detector systematic
  // holeice cached
  double cachedHoleIceConv;
  double cachedHoleIcePrompt;
  double cachedHoleIceAstro;

  // domeff cached
  double cachedDOMEffConv;
  double cachedDOMEffPrompt;
  double cachedDOMEffAstro;

  // year of the event
  unsigned int year;

  Event():
    //run(0),
    //event(0),
    //subevent(0),
    // LW/LI compatibility
    final_state_particle_0(LW::ParticleType::unknown),
    final_state_particle_1(LW::ParticleType::unknown),
    primaryType(LW::ParticleType::unknown),
    primaryEnergy(std::numeric_limits<float>::quiet_NaN()),
    primaryAzimuth(std::numeric_limits<float>::quiet_NaN()),
    primaryZenith(std::numeric_limits<float>::quiet_NaN()),
    totalColumnDepth(std::numeric_limits<float>::quiet_NaN()),
    intX(std::numeric_limits<float>::quiet_NaN()),
    intY(std::numeric_limits<float>::quiet_NaN()),
    oneWeight(0),
    number_of_generated_mc_events(std::numeric_limits<float>::quiet_NaN()),
    // meta-event quantities
    num_events(1),
    // reconstruction quantities
    energy(std::numeric_limits<float>::quiet_NaN()),
    zenith(std::numeric_limits<float>::quiet_NaN()),
    azimuth(std::numeric_limits<float>::quiet_NaN()),
    // which sample
    sample(golemfit::sampleTag::Sterile),
    // topology
    topology(0),
    // cached quantities
    cachedLivetime(0),
    cachedWeight(0),
    // flux caches
    cachedConvPionWeight(0),
    cachedConvKaonWeight(0),
    cachedConvWeight(0),
    cachedPromptWeight(0),
    cachedAstroMuWeight(0),
    // barr parameterers caches
    cachedBarrModWP(0),
    cachedBarrModWM(0),
    cachedBarrModYP(0),
    cachedBarrModYM(0),
    cachedBarrModZP(0),
    cachedBarrModZM(0),
    // holeice cache
    cachedHoleIceConv(0),
    cachedHoleIcePrompt(0),
    cachedHoleIceAstro(0),
    // dom efficiency cache
    cachedDOMEffConv(0),
    cachedDOMEffPrompt(0),
    cachedDOMEffAstro(0),
    year(0)
  {}

  babyEvent spitBabyEvent() const {
    assert(not std::isnan(oneWeight));
    babyEvent bevt;
    // LW/LI compatibility
    bevt.primaryType = primaryType;
    bevt.primaryEnergy = primaryEnergy;
    bevt.primaryAzimuth = primaryAzimuth;
    bevt.primaryZenith = primaryZenith;
    // reconstructed quantities
    bevt.energy = energy;
    bevt.zenith = zenith;
    bevt.azimuth = azimuth;
    // weight
    bevt.oneWeight = oneWeight;
    // topology
    bevt.topology = topology;
    return bevt;
  }

  bool operator==(const Event& e) {
    return  (((final_state_particle_0 == e.final_state_particle_0)) &&
            ((final_state_particle_1 == e.final_state_particle_1)) &&
            ((primaryType == e.primaryType)) &&
            ((primaryEnergy == e.primaryEnergy) || (std::isnan(primaryEnergy) && std::isnan(e.primaryEnergy))) &&
            ((primaryAzimuth == e.primaryAzimuth) || (std::isnan(primaryAzimuth) && std::isnan(e.primaryAzimuth))) &&
            ((primaryZenith == e.primaryZenith) || (std::isnan(primaryZenith) && std::isnan(e.primaryZenith))) &&
            ((totalColumnDepth == e.totalColumnDepth) || (std::isnan(totalColumnDepth) && std::isnan(e.totalColumnDepth))) &&
            ((intX == e.intX) || (std::isnan(intX) && std::isnan(e.intX))) &&
            ((intY == e.intY) || (std::isnan(intY) && std::isnan(e.intY))) &&
            ((oneWeight == e.oneWeight) || (std::isnan(oneWeight) && std::isnan(e.oneWeight))) &&
            ((number_of_generated_mc_events == e.number_of_generated_mc_events) || (std::isnan(number_of_generated_mc_events) && std::isnan(e.number_of_generated_mc_events))) &&
            ((num_events == e.num_events) || (std::isnan(num_events) && std::isnan(e.num_events))) &&
            ((energy == e.energy) || (std::isnan(energy) && std::isnan(e.energy))) &&
            ((zenith == e.zenith) || (std::isnan(zenith) && std::isnan(e.zenith))) &&
            ((azimuth == e.azimuth) || (std::isnan(azimuth) && std::isnan(e.azimuth))) &&
            ((sample == e.sample)) &&
            ((topology == e.topology) || (std::isnan(topology) && std::isnan(e.topology))) &&
            ((cachedLivetime == e.cachedLivetime) || (std::isnan(cachedLivetime) && std::isnan(e.cachedLivetime))) &&
            ((cachedConvPionWeight == e.cachedConvPionWeight) || (std::isnan(cachedConvPionWeight) && std::isnan(e.cachedConvPionWeight))) &&
            ((cachedConvKaonWeight == e.cachedConvKaonWeight) || (std::isnan(cachedConvKaonWeight) && std::isnan(e.cachedConvKaonWeight))) &&
            ((cachedConvWeight == e.cachedConvWeight) || (std::isnan(cachedConvWeight) && std::isnan(e.cachedConvWeight))) &&
            ((cachedPromptWeight == e.cachedPromptWeight) || (std::isnan(cachedPromptWeight) && std::isnan(e.cachedPromptWeight))) &&
            ((cachedAstroMuWeight == e.cachedAstroMuWeight) || (std::isnan(cachedAstroMuWeight) && std::isnan(e.cachedAstroMuWeight))) &&
            ((cachedWeight == e.cachedWeight) || (std::isnan(cachedWeight) && std::isnan(e.cachedWeight))) &&
            ((cachedBarrModWP == e.cachedBarrModWP) || (std::isnan(cachedBarrModWP) && std::isnan(e.cachedBarrModWP))) &&
            ((cachedBarrModWM == e.cachedBarrModWM) || (std::isnan(cachedBarrModWM) && std::isnan(e.cachedBarrModWM))) &&
            ((cachedBarrModYP == e.cachedBarrModYP) || (std::isnan(cachedBarrModYP) && std::isnan(e.cachedBarrModYP))) &&
            ((cachedBarrModYM == e.cachedBarrModYM) || (std::isnan(cachedBarrModYM) && std::isnan(e.cachedBarrModYM))) &&
            ((cachedBarrModZP == e.cachedBarrModZP) || (std::isnan(cachedBarrModZP) && std::isnan(e.cachedBarrModZP))) &&
            ((cachedBarrModZM == e.cachedBarrModZM) || (std::isnan(cachedBarrModZM) && std::isnan(e.cachedBarrModZM))) &&
            ((cachedHoleIceConv == e.cachedHoleIceConv) || (std::isnan(cachedHoleIceConv) && std::isnan(e.cachedHoleIceConv))) &&
            ((cachedHoleIcePrompt == e.cachedHoleIcePrompt) || (std::isnan(cachedHoleIcePrompt) && std::isnan(e.cachedHoleIcePrompt))) &&
            ((cachedHoleIceAstro== e.cachedHoleIceAstro) || (std::isnan(cachedHoleIceAstro) && std::isnan(e.cachedHoleIceAstro))) &&
            ((cachedDOMEffConv == e.cachedDOMEffConv) || (std::isnan(cachedDOMEffConv) && std::isnan(e.cachedDOMEffConv))) &&
            ((cachedDOMEffPrompt == e.cachedDOMEffPrompt) || (std::isnan(cachedDOMEffPrompt) && std::isnan(e.cachedDOMEffPrompt))) &&
            ((cachedDOMEffAstro == e.cachedDOMEffAstro) || (std::isnan(cachedDOMEffAstro) && std::isnan(e.cachedDOMEffAstro))) &&
            ((year == e.year) || (std::isnan(year) && std::isnan(e.year)))
            );
  }

  double Get_MetaWeight() const {
    //double meta_weight = oneWeight/number_of_generated_mc_events;
    double meta_weight = (cachedConvWeight + cachedPromptWeight + cachedAstroMuWeight);
    return meta_weight;
  }

  void operator+=(const Event& e) {
    // if the meta_weight is not 1, but its the oneWeight
    // the weighted fast mode is been turn on. This is a
    // super sharp knive. Appropiate magic needs to be done in the
    // division guy. CA.
    double meta_weight = e.Get_MetaWeight();
    assert((not std::isnan(meta_weight)) and (meta_weight >= 0.0));

    num_events += e.num_events;
    energy += e.energy*meta_weight;
    zenith += e.zenith*meta_weight;
    azimuth += e.azimuth*meta_weight;
    primaryEnergy += e.primaryEnergy*meta_weight;
    primaryZenith += e.primaryZenith*meta_weight;
    primaryAzimuth += e.primaryAzimuth*meta_weight;
    totalColumnDepth += e.totalColumnDepth*meta_weight;
    intX += e.intX*meta_weight;
    intY += e.intY*meta_weight;
    oneWeight += e.oneWeight;
    cachedLivetime += e.cachedLivetime*meta_weight;
    // we accumulate the rates
    cachedWeight += e.cachedWeight;
    cachedConvKaonWeight += e.cachedConvKaonWeight;
    cachedConvPionWeight += e.cachedConvPionWeight;
    cachedConvWeight += e.cachedConvWeight;
    cachedPromptWeight += e.cachedPromptWeight;
    cachedAstroMuWeight += e.cachedAstroMuWeight;
    cachedBarrModWP += e.cachedBarrModWP;
    cachedBarrModWM += e.cachedBarrModWM;
    cachedBarrModYP += e.cachedBarrModYP;
    cachedBarrModYM += e.cachedBarrModYM;
    cachedBarrModZP += e.cachedBarrModZP;
    cachedBarrModZM += e.cachedBarrModZM;
    // we average the expectation from systematics variation in reco bin
    cachedHoleIceConv = log10(pow(10.,cachedHoleIceConv) + pow(10.,e.cachedHoleIceConv)*meta_weight);
    cachedHoleIcePrompt = log10(pow(10.,cachedHoleIcePrompt) + pow(10.,e.cachedHoleIcePrompt)*meta_weight);
    cachedHoleIceAstro = log10(pow(10.,cachedHoleIceAstro) + pow(10.,e.cachedHoleIceAstro)*meta_weight);
    cachedDOMEffConv = log10(pow(10.,cachedDOMEffConv) + pow(10.,e.cachedDOMEffConv)*meta_weight);
    cachedDOMEffPrompt = log10(pow(10.,cachedDOMEffPrompt) + pow(10.,e.cachedDOMEffPrompt)*meta_weight);
    cachedDOMEffAstro = log10(pow(10.,cachedDOMEffAstro) + pow(10.,e.cachedDOMEffAstro)*meta_weight);
    year += e.year*meta_weight;
  }

  void operator/=(double d) {
    energy /= d;
    zenith /= d;
    azimuth /= d;
    primaryEnergy /= d;
    primaryZenith /= d;
    primaryAzimuth /= d;
    totalColumnDepth /= d;
    intX /= d;
    intY /= d;
    cachedLivetime /= d;
    year /= d;
    cachedHoleIceConv = log10(pow(10.,cachedHoleIceConv)/d);
    cachedHoleIcePrompt = log10(pow(10.,cachedHoleIcePrompt)/d);
    cachedHoleIceAstro = log10(pow(10.,cachedHoleIceAstro)/d);
    cachedDOMEffConv = log10(pow(10.,cachedDOMEffConv)/d);
    cachedDOMEffPrompt = log10(pow(10.,cachedDOMEffPrompt)/d);
    cachedDOMEffAstro= log10(pow(10.,cachedDOMEffAstro)/d);
  }

  bool check(bool check_, Level level) const{
    if (!check_) {
      return true;
    } else if ( level == Level::generation){
      return(!std::isnan(primaryEnergy) && !std::isnan(primaryZenith) && !std::isnan(primaryAzimuth)
         && !std::isnan(totalColumnDepth)
         && !std::isnan(intX) && !std::isnan(intY)
         && primaryType!=LW::ParticleType::unknown
         && final_state_particle_0!=LW::ParticleType::unknown && final_state_particle_0!=LW::ParticleType::unknown);
    } else {
      return false;
    }
  }
};

// input output
std::ostream& operator<<(std::ostream& os, const Event& e);

template<typename dataType>
Event combine_events(const std::vector<dataType>& events) {
  golemfit::fastmode::detail::dereference<dataType> dr;
  Event meta_e(events[0]);
  double w = meta_e.Get_MetaWeight();
  meta_e /= (1./meta_e.Get_MetaWeight());
  for(unsigned int i=1; i<events.size(); ++i) {
    //if(dr(events[i]).oneWeight>0){
      meta_e += events[i];
      w += dr(events[i]).Get_MetaWeight();
    //}
  }
  if(w < 0)
    throw std::runtime_error("Spot a negative weight while combining events. Kill.");
  if(w > 0){
    meta_e /= w;
  }
  //std::cout << w << std::endl;
  //std::cout << meta_e << std::endl;
  // note that if w = 0; none of the events have weight and we just return the guy.
  return meta_e;
}

using json=nlohmann::json;

namespace nlohmann {
    template<>
    struct adl_serializer<LW::ParticleType> {
        static void to_json(json& j, const LW::ParticleType& p) {
            j = int(p);
        }

        static void from_json(const json& j, LW::ParticleType& p) {
            p = LW::ParticleType(j.get<int>());
        }
    };
}

struct double_proxy {
    template<typename T>
    double operator()(const T& t) {
        if(t.is_null())
            return std::numeric_limits<double>::quiet_NaN();
        else
            return t.template get<double>();
    }
};

struct float_proxy {
    template<typename T>
    double operator()(const T& t) {
        if(t.is_null())
            return std::numeric_limits<float>::quiet_NaN();
        else
            return t.template get<float>();
    }
};

void to_json(json& j, const Event& e);
void from_json(const json& j, Event& e);

namespace H5Load {

namespace sterile {

template<typename CallbackType>
void readFile(const std::string& filePath, CallbackType action){
    using namespace phys_tools::cts;
    using namespace phys_tools::tableio;
    H5File h5file(filePath);
    if(!h5file)
        throw std::runtime_error("Unable to open "+filePath);
    std::set<std::string> tables;
    H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
    if(tables.empty())
        throw std::runtime_error(filePath+" contains no tables");
    std::map<RecordID,Event> intermediateData;

    using double_value_field = TableRow<field<double,CTS("value")>>;

    if(tables.count("PrimaryType")){
        readTable<double_value_field>(h5file,"PrimaryType", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryType=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    if(tables.count("FinalStateX")){
        readTable<double_value_field>(h5file,"FinalStateX", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.intX=p.get<CTS("value")>();
            });
    }

    if(tables.count("FinalStateY")){
        readTable<double_value_field>(h5file,"FinalStateY", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.intY=p.get<CTS("value")>();
            });
    }

    if(tables.count("FinalType0")){
        readTable<double_value_field>(h5file,"FinalType0", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.final_state_particle_0=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    if(tables.count("FinalType1")){
        readTable<double_value_field>(h5file,"FinalType1", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.final_state_particle_1=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    /*
    if(tables.count("ImpactParameter")){
        readTable<double_value_field>(h5file,"ImpactParameter", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.r=(p.get<CTS("value")>());
            });
    }
    */

    if(tables.count("MuExAzimuth")){
        readTable<double_value_field>(h5file,"MuExAzimuth", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.azimuth=(p.get<CTS("value")>());
            });
    }

    if(tables.count("MuExZenith")){
        readTable<double_value_field>(h5file,"MuExZenith", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.zenith=(p.get<CTS("value")>());
            });
    }

    if(tables.count("MuExEnergy")){
        readTable<double_value_field>(h5file,"MuExEnergy", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.energy=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuAzimuth")){
        readTable<double_value_field>(h5file,"NuAzimuth", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryAzimuth=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuZenith")){
        readTable<double_value_field>(h5file,"NuZenith", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryZenith=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuEnergy")){
        readTable<double_value_field>(h5file,"NuEnergy", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryEnergy=(p.get<CTS("value")>());
            });
    }

    if(tables.count("TotalColumnDepth")){
        readTable<double_value_field>(h5file,"TotalColumnDepth", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.totalColumnDepth=(p.get<CTS("value")>());
            });
    }

    for(std::map<RecordID,Event>::value_type& item : intermediateData)
        action(item.first,item.second);

    intermediateData.clear();
}

} // close io sterile namespace

} // close h5load input namespace

#endif //EVENT_
