#include "Event.h"
#include "json.hpp"

std::ostream& operator<<(std::ostream& os, const Event& e){
    std::cout << " energy: " << e.energy << '\n'
    << " zenith: " << e.zenith << '\n'
    << " azimuth: " << e.azimuth << '\n'
    << " topology: " << (int)e.topology << '\n'
    << " primaryType: " << (int)e.primaryType << '\n'
    << " primaryZenith: " << e.primaryZenith << '\n'
    << " primaryAzimuth: " << e.primaryAzimuth << '\n'
    << " primaryEnergy: " << e.primaryEnergy << '\n'
    << " cachedConvPionWeight: " << e.cachedConvPionWeight << '\n'
    << " cachedConvKaonWeight: " << e.cachedConvKaonWeight << '\n'
    << " cachedPromptWeight: " << e.cachedPromptWeight << '\n'
    << " cachedBarrModWP: " << e.cachedBarrModWP << '\n'
    << " cachedBarrModWM: " << e.cachedBarrModWM << '\n'
    << " cachedBarrModYP: " << e.cachedBarrModYP << '\n'
    << " cachedBarrModYM: " << e.cachedBarrModYM << '\n'
    << " cachedBarrModZP: " << e.cachedBarrModZP << '\n'
    << " cachedBarrModZM: " << e.cachedBarrModZM << '\n'
    << " cachedDOMEffConv: " << e.cachedDOMEffConv<< '\n'
    << " cachedDOMEffPrompt: " << e.cachedDOMEffPrompt << '\n'
    << " cachedDOMEffAstro: " << e.cachedDOMEffAstro << '\n'
    << " cachedHoleIceConv: " << e.cachedHoleIceConv << '\n'
    << " cachedHoleIcePrompt: " << e.cachedHoleIcePrompt << '\n'
    << " cachedHoleIceAstro: " << e.cachedHoleIceAstro << '\n'
    << " OneWeight-TimesLifetime: " << e.oneWeight<< '\n';
    return(os);
}

double differenceAngle(double zenith1, double azimuth1, double zenith2, double azimuth2){
	using boost::math::constants::pi;

	double cad=cos(azimuth1-azimuth2);
	double czd=cos(zenith1-zenith2);
	double czs=cos(zenith1+zenith2);
	double dot=(cad*(czd-czs)+czd+czs)/2;
	if(dot>1.)
		return(0.0);
	if(dot<-1.)
		return(pi<double>());
	return(acos(dot));
}

herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data){
	std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
	items->insert(member_name);
	return(0);
}

void to_json(json& j, const Event& e) {
    j = {
        {"final_state_particle_0", e.final_state_particle_0},
        {"final_state_particle_1", e.final_state_particle_1},
        {"primaryType", e.primaryType},
        {"primaryEnergy", e.primaryEnergy},
        {"primaryAzimuth", e.primaryAzimuth},
        {"primaryZenith", e.primaryZenith},
        {"totalColumnDepth", e.totalColumnDepth},
        {"intX", e.intX},
        {"intY", e.intY},
        {"oneWeight", e.oneWeight},
        {"number_of_generated_mc_events", e.number_of_generated_mc_events},
        {"num_events", e.num_events},
        {"energy", e.energy},
        {"zenith", e.zenith},
        {"azimuth", e.azimuth},
        {"sample", e.sample},
        {"topology", e.topology},
        {"cachedLivetime", e.cachedLivetime},
        {"cachedConvPionWeight", e.cachedConvPionWeight},
        {"cachedConvKaonWeight", e.cachedConvKaonWeight},
        {"cachedPromptWeight", e.cachedPromptWeight},
        {"cachedWeight", e.cachedWeight},
        {"cachedBarrModWP", e.cachedBarrModWP},
        {"cachedBarrModWM", e.cachedBarrModWM},
        {"cachedBarrModYP", e.cachedBarrModYP},
        {"cachedBarrModYM", e.cachedBarrModYM},
        {"cachedBarrModZP", e.cachedBarrModZP},
        {"cachedBarrModZM", e.cachedBarrModZM},
        {"cachedHoleIceConv", e.cachedHoleIceConv},
        {"cachedHoleIcePrompt", e.cachedHoleIcePrompt},
        {"cachedDOMEffConv", e.cachedDOMEffConv},
        {"cachedDOMEffPrompt", e.cachedDOMEffPrompt}
    };
}

void from_json(const json& j, Event& e) {
    double_proxy dp;
    float_proxy fp;
    e.final_state_particle_0 = j.at("final_state_particle_0").get<std::remove_reference<decltype(e.final_state_particle_0)>::type>();
    e.final_state_particle_1 = j.at("final_state_particle_1").get<std::remove_reference<decltype(e.final_state_particle_1)>::type>();
    e.primaryType = j.at("primaryType").get<std::remove_reference<decltype(e.primaryType)>::type>();
    e.primaryEnergy = dp(j.at("primaryEnergy"));
    e.primaryAzimuth = dp(j.at("primaryAzimuth"));
    e.primaryZenith = dp(j.at("primaryZenith"));
    e.totalColumnDepth = dp(j.at("totalColumnDepth"));
    e.intX = dp(j.at("intX"));
    e.intY = dp(j.at("intY"));
    e.oneWeight = fp(j.at("oneWeight"));
    e.number_of_generated_mc_events = fp(j.at("number_of_generated_mc_events"));
    e.num_events = j.at("num_events").get<std::remove_reference<decltype(e.num_events)>::type>();
    e.energy = fp(j.at("energy"));
    e.zenith = fp(j.at("zenith"));
    e.azimuth = fp(j.at("azimuth"));
    e.sample = j.at("sample").get<std::remove_reference<decltype(e.sample)>::type>();
    e.topology = j.at("topology").get<std::remove_reference<decltype(e.topology)>::type>();
    e.cachedLivetime = j.at("cachedLivetime").get<std::remove_reference<decltype(e.cachedLivetime)>::type>();
    e.cachedConvPionWeight = j.at("cachedConvPionWeight").get<std::remove_reference<decltype(e.cachedConvPionWeight)>::type>();
    e.cachedConvKaonWeight = j.at("cachedConvKaonWeight").get<std::remove_reference<decltype(e.cachedConvKaonWeight)>::type>();
    e.cachedPromptWeight = j.at("cachedPromptWeight").get<std::remove_reference<decltype(e.cachedPromptWeight)>::type>();
    e.cachedWeight = j.at("cachedWeight").get<std::remove_reference<decltype(e.cachedWeight)>::type>();

    e.cachedBarrModWP = j.at("cachedBarrModWP").get<std::remove_reference<decltype(e.cachedBarrModWP)>::type>();
    e.cachedBarrModWM = j.at("cachedBarrModWM").get<std::remove_reference<decltype(e.cachedBarrModWM)>::type>();
    e.cachedBarrModYP = j.at("cachedBarrModYP").get<std::remove_reference<decltype(e.cachedBarrModYP)>::type>();
    e.cachedBarrModYM = j.at("cachedBarrModYM").get<std::remove_reference<decltype(e.cachedBarrModYM)>::type>();
    e.cachedBarrModZP = j.at("cachedBarrModZP").get<std::remove_reference<decltype(e.cachedBarrModZP)>::type>();
    e.cachedBarrModZM = j.at("cachedBarrModZM").get<std::remove_reference<decltype(e.cachedBarrModZM)>::type>();

    e.cachedHoleIceConv = j.at("cachedHoleIceConv").get<std::remove_reference<decltype(e.cachedHoleIceConv)>::type>();
    e.cachedHoleIcePrompt = j.at("cachedHoleIcePrompt").get<std::remove_reference<decltype(e.cachedHoleIcePrompt)>::type>();

    e.cachedDOMEffConv = j.at("cachedDOMEffConv").get<std::remove_reference<decltype(e.cachedDOMEffConv)>::type>();
    e.cachedDOMEffPrompt = j.at("cachedDOMEffPrompt").get<std::remove_reference<decltype(e.cachedDOMEffPrompt)>::type>();
}
