#ifndef RUN_SPEC
#define RUN_SPEC

#include <LeptonWeighter/Event.h>
#include <LeptonWeighter/ParticleType.h>
#include <LeptonWeighter/Generator.h>

#include "GolemEnumDefinitions.h"

namespace golemfit {

struct MCSet {
  const MCType type;
	const std::string path;
  const std::string filename;
  const std::pair<bool, unsigned int> split;
  const double number_of_generated_events;
  const double unshadowedFraction ;// aka domefficiency
  const double holeiceForward;// aka p2
  const double anisotropyScale;
  std::vector<std::shared_ptr<LW::Generator>> generators;
  // constructor
	MCSet(MCType type, const std::string& p,const std::string& filename, const std::pair<bool, unsigned int> split, double number_of_generated_events, double unshadowedFraction, double holeiceForward, double anisotropyScale, std::vector<std::shared_ptr<LW::Generator>> g):
	type(type),path(p),filename(filename),split(split),number_of_generated_events(number_of_generated_events),unshadowedFraction(unshadowedFraction),holeiceForward(holeiceForward),anisotropyScale(anisotropyScale),generators(g){}
  // copy constructor
	MCSet(MCSet& other):type(other.type),path(other.path),
  filename(other.filename),split(other.split),number_of_generated_events(other.number_of_generated_events),
  unshadowedFraction(other.unshadowedFraction),holeiceForward(other.holeiceForward),anisotropyScale(other.anisotropyScale),generators(other.generators){}
  // copy constructor
	MCSet(const MCSet& other):type(other.type),path(other.path),
  filename(other.filename),split(other.split),number_of_generated_events(other.number_of_generated_events),
  unshadowedFraction(other.unshadowedFraction),holeiceForward(other.holeiceForward),anisotropyScale(other.anisotropyScale),generators(other.generators){}
  // move constructor
	MCSet(MCSet&& other):type(std::move(other.type)),path(std::move(other.path)),
  filename(std::move(other.filename)),split(std::move(other.split)),number_of_generated_events(std::move(other.number_of_generated_events)),
  unshadowedFraction(std::move(other.unshadowedFraction)),holeiceForward(std::move(other.holeiceForward)),anisotropyScale(std::move(other.anisotropyScale)),generators(std::move(other.generators)){}
};

} // close golemfit namespace

#endif
