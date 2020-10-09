#ifndef ANALYSISWEIGHTING_H
#define ANALYSISWEIGHTING_H

#include <cassert>
#include <type_traits>
#include <limits>

//#include <LeptonWeighter/particleType.h>
//#include "../Likelihood/autodiff.h"
//#include "weighting.h"
#include <PhysTools/likelihood/weighting.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/autodiff.h>
#include <NewNuFlux/NewNuFlux.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include <nuFATE/nuFATE.h>
#include <math.h>

#include "GolemParameters.h"
#include "GolemTools.h"
#include "utils.h"

//================================================================================
// FUNDAMENTAL WEIGHTERS
//================================================================================

template<typename T, typename Event, typename U>
struct cachedValueWeighter : public phys_tools::GenericWeighter<cachedValueWeighter<T,Event,U>>{
    private:
        U Event::* cachedPtr;
    public:
        using result_type=T;
        cachedValueWeighter(U Event::* ptr):cachedPtr(ptr){}
        result_type operator()(const Event& e) const{
            return(result_type(e.*cachedPtr));
        }
};

template<typename T, typename Event, typename U, typename K>
struct cachedMapValueWeighter : public phys_tools::GenericWeighter<cachedMapValueWeighter<T,Event,U,K>>{
    private:
        U Event::* cachedPtr;
        const K key;
    public:
        using result_type=T;
        cachedMapValueWeighter(U Event::* ptr,K key):cachedPtr(ptr),key(key){}
        result_type operator()(const Event& e) const{
            return(result_type((e.*cachedPtr).at(key)));
        }
};

template<typename T, typename Event>
struct FunctionWeighter : public phys_tools::GenericWeighter<FunctionWeighter<T,Event>>{
    private:
        std::function<T(const Event&)> func;
    public:
        using result_type=T;
        FunctionWeighter(std::function<T(const Event&)> f):func(f){}
        result_type operator()(const Event& e) const{
            return(result_type(func(e)));
        }
};

template<typename T, typename E>
FunctionWeighter<T,E> makeFunctionWeighter(std::function<T(const E&)> f){
    return(FunctionWeighter<T,E>(f));
}

template<typename Event, typename T>
struct boxEnergyFilter : public phys_tools::GenericWeighter<boxEnergyFilter<Event,T>>{
    private:
        double min, max;
    public:
        using result_type=T;
        boxEnergyFilter(double min, double max):min(min),max(max){}

        result_type operator()(const Event& e) const{
            return((e.primaryEnergy>min && e.primaryEnergy<max) ? 1 : 0);
        }
};

//================================================================================
// ADHOC TABLE WEIGHTER WEIGHTERS
//================================================================================

template<unsigned int HistogramDim>
struct tableWeighter: public phys_tools::GenericWeighter<tableWeighter<HistogramDim>>{
    private:
        unsigned int get_meta_index(const std::vector<unsigned int>& indices, const std::vector<unsigned int>& bin_per_dimension) const {
            assert(indices.size()==bin_per_dimension.size());
            unsigned int index = 0;
            for(unsigned i = 0; i < indices.size(); i++){
                unsigned int meta_size = 1;
                for(unsigned int j = i; j < indices.size(); j++)
                    meta_size *= bin_per_dimension[j];
                index += meta_size*indices[i];
            }
            return index;
        }
        std::vector<unsigned int> get_indices(std::vector<double> coordinates) const {
          std::vector<unsigned int> indices(HistogramDim);
          for(unsigned int i=0; i< HistogramDim; i++){
            auto it=std::lower_bound(histogram_edges_array[i].begin(),histogram_edges_array[i].end(),coordinates[i]);
            if(it==histogram_edges_array[i].end())
              throw std::runtime_error("This code does not extrapolate.");
            if(it!=histogram_edges_array[i].begin())
              it--;
            indices[i] = std::distance(histogram_edges_array[i].begin(),it);
          }
          return(indices);
        }
    protected:
        std::vector<double> eta_pivots;
        std::vector<nusquids::marray<double,HistogramDim>> table;
        const std::vector<std::vector<double>> histogram_edges_array;
    public:
        tableWeighter(std::vector<double> eta_pivots, std::vector<nusquids::marray<double,HistogramDim>> table, std::vector<std::vector<double>> histogram_edges_array):
          eta_pivots(eta_pivots),table(table),histogram_edges_array(histogram_edges_array)
        {
            assert(eta_pivots.size() == table.size());
            assert(HistogramDim >1);
            assert(histogram_edges_array.size() == HistogramDim);
        }
        tableWeighter(std::vector<double> eta_pivots, std::vector<std::string> list_of_filenames, std::vector<std::vector<double>> histogram_edges_array):
          eta_pivots(eta_pivots),histogram_edges_array(histogram_edges_array)
        {
            assert(eta_pivots.size() == list_of_filenames.size());
            assert(HistogramDim >1);
            assert(histogram_edges_array.size() == HistogramDim);

            unsigned int expected_file_length = 1;
            std::vector<size_t> bins_per_dimension;
            for(unsigned int i = 0; i < histogram_edges_array.size(); i++){
                expected_file_length *= (histogram_edges_array[i].size()-1); // number of bins per dimension
                bins_per_dimension.push_back((histogram_edges_array[i].size()-1));
            }

            if(expected_file_length == 0 or expected_file_length == 1)
                throw std::runtime_error("Error::tableWeighter::tableWeighter: expected histogram_edges_array dimension do not make sense");

            for(std::string filename : list_of_filenames){
                auto data_guy = nusquids::quickread(filename);
                if(data_guy.extent(0) != expected_file_length)
                    throw std::runtime_error("Error::tableWeighter::tableWeighter: file " + filename + " has " + std::to_string(data_guy.extent(0)) + " but it should be " + std::to_string(expected_file_length));

                nusquids::marray<double,HistogramDim> tabela {bins_per_dimension};
                std::vector<unsigned int> indices {HistogramDim};
                for(unsigned int j = 0; j < indices.size(); j++){
                    for(; indices[j] < bins_per_dimension[j]; indices[j]++){
                        unsigned int meta_index = get_meta_index(indices,bins_per_dimension);
                        tabela[indices] = data_guy[meta_index][HistogramDim];
                    }
                }
                table.push_back(tabela);
            }
        }

        nusquids::marray<double,2> GetHistogramEdges() const {
          return((nusquids::marray<double,2> const)histogram_edges_array);
        }

        std::vector<nusquids::marray<double,HistogramDim>> GetTable() const {
          return(table);
        }

        std::vector<double> GetParameterPivots() const {
          return(eta_pivots);
        }

        double GetInterpolatedValue(std::vector<double> coordinates, double eta) const {
          if(coordinates.size() != HistogramDim)
            throw std::runtime_error("Incorrect number of coordinates supplied.");
          std::vector<unsigned int> indices = get_indices(coordinates);;
          return GetInterpolatedValue(indices,eta);
        }

        double GetInterpolatedValue(std::vector<unsigned int> indices, double eta) const {
          if(indices.size() != HistogramDim)
            throw std::runtime_error("Incorrect number of indices supplied.");
          auto it=std::lower_bound(eta_pivots.begin(),eta_pivots.end(),eta);
          if(it==eta_pivots.end())
            throw std::runtime_error("This code does not extrapolate.");
          if(it!=eta_pivots.begin())
            it--;
          size_t eta_index=std::distance(eta_pivots.begin(),it);

          auto & lower_array = table[eta_index];
          auto & higher_array = table[eta_index+1];

          return lower_array[indices] + (higher_array[indices]-lower_array[indices])*(eta-eta_pivots[eta_index])/(eta_pivots[eta_index+1] - eta_pivots[eta_index]);
        }

        double GetGradient(std::vector<double> coordinates, double eta) const {
          if(coordinates.size() != HistogramDim)
            throw std::runtime_error("Incorrect number of coordinates supplied.");
          std::vector<unsigned int> indices = get_indices(coordinates);;
          return GetGradient(indices,eta);
        }

        double GetGradient(std::vector<unsigned int> indices, double eta) const {
          if(indices.size() != HistogramDim)
            throw std::runtime_error("Incorrect number of indices supplied.");
          auto it=std::lower_bound(eta_pivots.begin(),eta_pivots.end(),eta);
          if(it==eta_pivots.end())
            throw std::runtime_error("This code does not extrapolate.");
          if(it!=eta_pivots.begin())
            it--;
          size_t eta_index=std::distance(eta_pivots.begin(),it);

          auto & lower_array = table[eta_index];
          auto & higher_array = table[eta_index+1];

          return (higher_array[indices]-lower_array[indices])/(eta_pivots[eta_index+1] - eta_pivots[eta_index]);
        }
};

template<typename Event, typename DataType>
struct holeIceTableWeighter: public tableWeighter<2>{
    private:
      const DataType eta;
    public:
      using table_struct=std::vector<nusquids::marray<double,2>>;
      holeIceTableWeighter(std::vector<double> eta_pivot, table_struct& table, std::vector<std::vector<double>>& histogram_edges_array, DataType eta):
        tableWeighter<2>(eta_pivot, table, histogram_edges_array), eta(eta){}

      using result_type=double;
      result_type operator()(const Event& e) const{
          return(GetInterpolatedValue({log10(e.energy),cos(e.zenith)},eta))/e.cachedHoleIceness;
      }
};

template<typename Event, unsigned int DerivativeDim>
struct holeIceTableWeighter<Event,phys_tools::autodiff::FD<DerivativeDim>>: public tableWeighter<2>{
     using result_type=phys_tools::autodiff::FD<DerivativeDim>;
    private:
      const result_type eta;
      unsigned int didx;
    public:
      using table_struct=std::vector<nusquids::marray<double,2>>;
      holeIceTableWeighter(std::vector<double> eta_pivot, table_struct& table, std::vector<std::vector<double>>& histogram_edges_array, result_type eta):
        tableWeighter<2>(eta_pivot,table, histogram_edges_array), eta(eta)
      {
        for(int i=0; i<DerivativeDim; i++){ //this will break if Dim is Dynamic
            if(eta.derivative(i)!=0){
                didx=i;
                break;
            }
        }
      }

      result_type operator()(const Event& e) const{
          double eta_ = eta.value();
          double cache=e.cachedHoleIceness;
          double rate_change = GetInterpolatedValue({log10(e.energy),cos(e.zenith)},eta_)/cache;
          double derivative = GetGradient({log10(e.energy),cos(e.zenith)},eta_)/cache;
          derivative*=eta.derivative(didx);
          result_type r{rate_change};
          r.setDerivative(derivative,didx);
          return(r);
      }
};

//================================================================================
// ATMOSPHERIC FLUX WEIGHTERS
//================================================================================

struct atmosNeutrinoFluxWeighter : public phys_tools::GenericWeighter<atmosNeutrinoFluxWeighter>{
    private:
        std::shared_ptr<NewNuFlux::FluxFunction> flux;
    public:
        atmosNeutrinoFluxWeighter(NewNuFlux::FluxFunction* f):flux(f){}
        atmosNeutrinoFluxWeighter(std::shared_ptr<NewNuFlux::FluxFunction> f):flux(f){}
        atmosNeutrinoFluxWeighter(boost::shared_ptr<NewNuFlux::FluxFunction> f):flux(to_std_ptr(f)){}

        using result_type=double;

        template<typename Event>
            result_type operator()(const Event& e) const{
                //factor of 2 to convert to per-flavor convention
                return(2*flux->getFlux((I3Particle::ParticleType)e.primaryType, e.primaryEnergy, cos(e.primaryZenith)));
            }

        std::shared_ptr<NewNuFlux::FluxFunction> get(){
            return(flux);
        }
};


namespace {
    template<typename Event, typename DataType>
        struct neuaneu : public phys_tools::GenericWeighter<neuaneu<Event,DataType>> {
            private:
                DataType ratio;
            public:
                neuaneu(DataType ratio): ratio(ratio) {};

                using result_type=DataType;
                result_type operator()(const Event& e) const{
                    //std::cout << (int) e.primaryType << std::endl;
                    //std::cout << (int) particleType::MuMinus << std::endl;
                    if(e.primaryType == LW::ParticleType::MuMinus or
                            e.primaryType == LW::ParticleType::EMinus  or
                            e.primaryType == LW::ParticleType::TauMinus )
                        return ratio;
                    else
                        return 2.-ratio;
                }
        };
}

//Weight particles and antiparticles differently
template<typename Event, typename T>
struct antiparticleWeighter : public phys_tools::GenericWeighter<antiparticleWeighter<Event,T>>{
    private:
        //1 = equal weighting of particles and antiparticles,
        //0 = zero weight for particles, double weight for antiparticles
        //2 = double weight for particles, zero weight for antiparticles
        T balance;
    public:
        using result_type=T;
        antiparticleWeighter(T b):
            balance(b){}

        result_type operator()(const Event& e) const{
            return((int)e.primaryType<0 ? balance : 2-balance);
        }
};


template<typename Event, typename DataType>
struct AtmosphericZenithVariationCorrectionFactor: public phys_tools::GenericWeighter<AtmosphericZenithVariationCorrectionFactor<Event,DataType>> {
    private:
        DataType delta;
        double c0;
        double c1;
        double e0;
        double e1;
    public:
        AtmosphericZenithVariationCorrectionFactor(DataType delta): delta(delta),c0(0.4),c1(100.),e0(369.),e1(11279.) {};

        using result_type=DataType;
        result_type operator()(const Event& e) const{
            double ct = cos(e.zenith);
            return 1.0 + (ct+c0)*delta*(1.+(e.energy-e0)/e1)/(1.+exp(-2.*c1*(ct+c0)));
        }
};


//Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
template<typename Event, typename T>
struct powerlawTiltWeighter : public phys_tools::GenericWeighter<powerlawTiltWeighter<Event,T>>{
    private:
        double medianEnergy;
        T deltaIndex;
    public:
        using result_type=T;
        powerlawTiltWeighter(double me, T dg):
            medianEnergy(me),deltaIndex(dg){}

        result_type operator()(const Event& e) const{
            result_type weight=pow(e.primaryEnergy/medianEnergy,-deltaIndex);
            return(weight);
        }
};

//================================================================================
// ASTROPHYSICAL FLUX WEIGHTERS
//================================================================================

struct piecewisePowerlawFlux{
    private:
        struct powerlaw{
            double eMin, eMax;
            double norm, power;

            bool containsEnergy(double e) const{ return(e>=eMin && e<eMax); }
            double operator()(double e) const{ return(norm*pow(e,power)); }
        };
        friend bool operator<(const powerlaw&, const powerlaw&);

        static bool overlap(const powerlaw& p1, const powerlaw& p2){
            if(p1<p2)
                return(p1.eMax<=p2.eMin);
            return(p2.eMax<=p1.eMin);
        }

        static bool below_max(double e, const powerlaw& p){ return(e<p.eMax); }

        void read(std::istream& is);

        std::vector<powerlaw> pieces;

    public:
        explicit piecewisePowerlawFlux(std::string path){
            std::ifstream file(path.c_str());
            if(!file)
                throw std::runtime_error("Couldn't open "+path+" for reading");
            read(file);
        }

        double getFlux(double e) const{
            auto it=std::upper_bound(pieces.begin(),pieces.end(),e,below_max);
            if(it==pieces.end())
                return(0); //treat flux as zero where it is undefined
            if(!it->containsEnergy(e))
                return(0); //treat flux as zero where it is undefined
            return((*it)(e));
        }
};

template<typename T>
struct powerlawWeighter : public phys_tools::GenericWeighter<powerlawWeighter<T>>{
    private:
        T index;
    public:
        using result_type=T;
        powerlawWeighter(T i):index(i){}

        template<typename Event>
            result_type operator()(const Event& e) const{
                return(pow((double)e.primaryEnergy,index));
            }
};

//Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
template<typename Event, typename T>
struct brokenPowerlawTiltWeighter : public phys_tools::GenericWeighter<brokenPowerlawTiltWeighter<Event,T>>{
    private:
        double medianEnergy;
        T deltaIndex_1;
        T deltaIndex_2;
        T sigma;

    public:
        using result_type=T;
        brokenPowerlawTiltWeighter(double me, T dg1, T dg2, T s):
            medianEnergy(me),deltaIndex_1(dg1),deltaIndex_2(dg2),sigma(s){}

        result_type operator()(const Event& e) const{
            //                                          only attempt to compute this if we have a primary energy!
            result_type weight=(e.primaryEnergy?pow(e.primaryEnergy/medianEnergy,-deltaIndex_1)+(1-sigma)*pow(e.primaryEnergy/medianEnergy,-deltaIndex_1-deltaIndex_2):1);
            return(weight);
        }
};

//Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
template<typename Event, typename T>
struct logParaboloidWeighter : public phys_tools::GenericWeighter<logParaboloidWeighter<Event,T>>{
    private:
        double breakEnergy;
        T alpha;
        T beta;
    public:
        using result_type=T;
        logParaboloidWeighter(double breakEnergy, T alpha, T beta):
            breakEnergy(breakEnergy),alpha(alpha),beta(beta){}

        result_type operator()(const Event& e) const{
            //                                          only attempt to compute this if we have a primary energy!
            result_type weight=pow(e.primaryEnergy/breakEnergy,-alpha -beta*log10(e.primaryEnergy/breakEnergy));
            return(weight);
        }
};

//Apply an exponential cutoff with a given characteristic energy
template<typename Event, typename T>
struct expCutoffWeighter : public phys_tools::GenericWeighter<expCutoffWeighter<Event,T>>{
    private:
        T cutoffEnergy;
    public:
        using result_type=T;
        expCutoffWeighter(T ce):
            cutoffEnergy(ce){}

        result_type operator()(const Event& e) const{
            //only attempt to compute this is we have a primary energy!
            result_type weight=(e.primaryEnergy?exp(-e.primaryEnergy/cutoffEnergy):1);
            return(weight);
        }
};

//================================================================================
// DOM EFFICIENCY WEIGHTER
//================================================================================

/*
namespace DOMEff3{
    using namespace photospline;
    using namespace phys_tools::autodiff;
    // Classic ChrisWeaver/Sterili-hack DOMefficiency variants. Start here. CA.
    using DOMMapType=std::map<unsigned int,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct DOMEffWeighter : public phys_tools::GenericWeighter<DOMEffWeighter<Event,DataType>>{
            private:
                DOMMapType domEffMap;
                DataType logEff;
            public:
                DOMEffWeighter(DOMMapType domEffMap, DataType deltaEff):
                    domEffMap(domEffMap),logEff(log10(0.9*(1.0+deltaEff))){}

                using result_type=double;
                result_type operator()(const Event& e) const{
                    float cache=e.cachedDOMEff;
                    double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};

                    auto domcorrection = domEffMap.find(e.year);
                    if( domcorrection == domEffMap.end() )
                        throw std::runtime_error("Dom efficiency correction for year "+std::to_string(e.year) + " not found.");
                    double rate=(*(*domcorrection).second)(coordinates);
                    return(pow(10.,rate-cache));
                }
        };

    template<typename Event,unsigned int Dim>
        struct DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                DOMMapType domEffMap;
                phys_tools::autodiff::FD<Dim> logEff;
                unsigned int didx;
            public:
                DOMEffWeighter(DOMMapType domEffMap, phys_tools::autodiff::FD<Dim> deltaEff):
                    domEffMap(domEffMap),
                    logEff(log10(0.9*(1.0+deltaEff)))
            {
                const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(deltaEff);
                for(int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(deltaEff.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                }
            }
                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                    float cache=e.cachedDOMEff;
                    double rate, derivative;
                    double coordinates[3]={log10(e.energy),cos(e.zenith),logEff.value()};
                    int centers[3];

                    auto domcorrection = domEffMap.find(e.year);
                    if( domcorrection == domEffMap.end() )
                        throw std::runtime_error("Dom efficiency correction for year "+std::to_string(e.year) + " not found.");
                    rate=(*(*domcorrection).second)(coordinates);
                    if(not (*(*domcorrection).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space");
                    derivative=(*(*domcorrection).second).ndsplineeval(coordinates,centers,2);
                    derivative*=logEff.derivative(didx);

                    result_type r(rate);
                    r.setDerivative(derivative,didx);
                    return(pow(10.,r-cache));
                }
        };

    //used for initializing the per-event dom-eff related caches
    template<typename Event>
        struct simpleEffRate{
            splinetable<>* rate2011;
            double logEff;
            simpleEffRate(splinetable<>* r2011, double eff):
                rate2011(r2011),logEff(log10(eff)){}

            void setCache(Event& e) const{
                double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};
                splinetable<>* rateTable=nullptr;
                rateTable = rate2011;
                   //switch(e.year){
                   //case 2011: rateTable=rate2011; break;
                   //default: assert(false && "Unexpected year");
                   //}
                e.cachedDOMEff=(*rateTable)(coordinates);
            }
        };

} // close namespace DOMEff3

namespace DOMEff4{
    // This namespace is for Chris new DOMEff variation implementation
    // where he moves the object to the event class.
    // We keep the previous implementation for simplifity. CA.
    template<typename Event, typename DataType>
        struct DOMEffWeighter : public phys_tools::GenericWeighter<DOMEffWeighter<Event,DataType>>{
            private:
                DataType eff;
                float Event::* cachedFactor;
                photospline::splinetable<>::evaluator* Event::* corrector;
            public:
                using result_type=DataType;
                DOMEffWeighter(DataType e, photospline::splinetable<>::evaluator* Event::* c, float Event::* f):
                    eff(e),cachedFactor(f),corrector(c){}
                result_type operator()(const Event& e) const{
#ifdef ALLOW_INVALID_DOMEFF
                    if(!(e.*corrector))
                        return(1.);
#endif
                    double coords[2]={log10(e.energy),eff};
                    double factor=(*(e.*corrector))(coords);
                    return(pow(10.,factor-(e.*cachedFactor)));
                }
        };

    template<typename Event, unsigned int Dim>
        struct DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                phys_tools::autodiff::FD<Dim> eff;
                float Event::* cachedFactor;
                photospline::splinetable<>::evaluator* Event::* corrector;
                unsigned int didx;
            public:
                using result_type=phys_tools::autodiff::FD<Dim>;
                DOMEffWeighter(result_type e, photospline::splinetable<>::evaluator* Event::* c, float Event::* f):
                    eff(e),cachedFactor(f),corrector(c){
                        const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(e);
                        for(int i=0; i<n; i++){
                            if(e.derivative(i)!=0){
                                didx=i;
                                break;
                            }
                        }
                    }
                result_type operator()(const Event& e) const{
#ifdef ALLOW_INVALID_DOMEFF
                    if(!(e.*corrector))
                        return(1.);
#endif
                    double coords[2]={log10(e.energy),eff.value()};
                    int centers[2];
                    (e.*corrector)->searchcenters(coords,centers);
                    double factor=(e.*corrector)->ndsplineeval(coords,centers,0);
                    double derivative=(e.*corrector)->ndsplineeval(coords,centers,1<<1);
                    result_type r(factor);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r-(e.*cachedFactor)));
                }
        };
} // close DOMEff4 namespace
*/

namespace DOMEff5{
    using namespace photospline;
    using namespace phys_tools::autodiff;
    using DOMMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct DOMEffWeighter : public phys_tools::GenericWeighter<DOMEffWeighter<Event,DataType>>{
            private:
                const DOMMapType& dom_efficiency_map_;
                const DataType dom_efficiency_;
                const golemfit::FluxComponent flux_component_;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const DOMMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                DOMEffWeighter(const DOMMapType& dom_efficiency_map_, DataType dom_efficiency_, golemfit::FluxComponent flux_component_):
                    dom_efficiency_map_(dom_efficiency_map_),
                    dom_efficiency_(dom_efficiency_),
                    flux_component_(flux_component_),
                    enforce_needed((not dom_efficiency_map_.empty()) and IsComponentInMap(dom_efficiency_map_,flux_component_))
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                    if(dom_efficiency_map_.empty() or not enforce_needed){
                      return result_type(1.);
                    }
                    double cache;
                    if(flux_component_ == golemfit::FluxComponent::atmConv)
                      cache = e.cachedDOMEffConv;
                    else if (flux_component_ == golemfit::FluxComponent::atmPrompt)
                      cache = e.cachedDOMEffPrompt;
                    else if (flux_component_ == golemfit::FluxComponent::diffuseAstro_mu)
                      cache = e.cachedDOMEffAstro;
                    else
                      throw std::runtime_error("DOM efficiency correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");

                    double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_};

                    auto domefficiencycorrection = dom_efficiency_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(domefficiencycorrection == dom_efficiency_map_.end())
                        throw std::runtime_error("DOM efficiency correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    double rate=(*(*domefficiencycorrection).second)(coordinates);

                    if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                      // this is a cow, ignore the event
                      return 0.0;
                    }

                    return(pow(10.,rate-cache));
                }
        };

    template<typename Event,unsigned int Dim>
        struct DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const DOMMapType& dom_efficiency_map_;
                const phys_tools::autodiff::FD<Dim> dom_efficiency_;
                const golemfit::FluxComponent flux_component_;
                unsigned int didx;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const DOMMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                DOMEffWeighter(const DOMMapType& dom_efficiency_map_, phys_tools::autodiff::FD<Dim> dom_efficiency_, golemfit::FluxComponent flux_component_):
                    dom_efficiency_map_(dom_efficiency_map_),
                    dom_efficiency_(dom_efficiency_),
                    flux_component_(flux_component_),
                    enforce_needed((not dom_efficiency_map_.empty()) and IsComponentInMap(dom_efficiency_map_,flux_component_))
                {
                    //const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(dom_efficiency_);
                    for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                        if(dom_efficiency_.derivative(i)!=0){
                            didx=i;
                            break;
                        }
                    }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                    if(dom_efficiency_map_.empty() or not enforce_needed){
                      return result_type(1.);
                    }

                    double cache;
                    if(flux_component_ == golemfit::FluxComponent::atmConv)
                      cache = e.cachedDOMEffConv;
                    else if (flux_component_ == golemfit::FluxComponent::atmPrompt)
                      cache = e.cachedDOMEffPrompt;
                    else if (flux_component_ == golemfit::FluxComponent::diffuseAstro_mu)
                      cache = e.cachedDOMEffAstro;
                    else
                      throw std::runtime_error("DOM efficiency correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");

                    double rate, derivative;
                    double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_.value()};
                    int centers[3];

                    auto domefficiencycorrection = dom_efficiency_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(domefficiencycorrection == dom_efficiency_map_.end())
                        throw std::runtime_error("DOM efficiency correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    rate=(*(*domefficiencycorrection).second)(coordinates);

                    if(not (*(*domefficiencycorrection).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: domeff spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(dom_efficiency_.value()) + ").");
                    derivative=(*(*domefficiencycorrection).second).ndsplineeval(coordinates,centers,1<<2);
                    derivative*=dom_efficiency_.derivative(didx);

                    if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                      // this is a cow, ignore the event
                      result_type r(0.0);
                      r.setDerivative(didx,0.0);
                      return r;
                    }

                    result_type r(rate);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r-cache));
                }
        };

    //used for initializing the per-event dom-eff related caches
    template<typename Event>
        struct DOMEfficiencySetter{
            const DOMMapType& dom_efficiency_map_;
            const double dom_efficiency_;
            DOMEfficiencySetter(const DOMMapType& dom_efficiency_map_, double dom_efficiency_):
                dom_efficiency_map_(dom_efficiency_map_),dom_efficiency_(dom_efficiency_){}

            void setCache(Event& e) const{
                assert(not dom_efficiency_map_.empty());
                double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_};
                for(auto element : dom_efficiency_map_){
                  if(element.first.second == static_cast<golemfit::Topology>(e.topology)){
                    double cache = (*element.second)(coordinates);
                    assert(std::isfinite(cache));
                    if(cache < -1.e5){
                      std::cout << e << std::endl;
                      throw std::runtime_error("Cows do not pass.");
                    }
                    if(cache == 0){
                      // the cache if the log of the rate; this cases the event to be ignored.
                      cache = -std::numeric_limits<float>::max();
                    }
                    if(element.first.first == golemfit::FluxComponent::atmConv)
                      e.cachedDOMEffConv = cache;
                    else if(element.first.first == golemfit::FluxComponent::atmPrompt)
                      e.cachedDOMEffPrompt = cache;
                    else if(element.first.first == golemfit::FluxComponent::diffuseAstro_mu)
                      e.cachedDOMEffAstro= cache;
                    else
                      throw std::runtime_error("Unimplemented DOMEff component");
                  }
                }
            }
        };
} // close DOMEff5 namespace

//================================================================================
// HOLEICE WEIGHTER
//================================================================================

namespace HoleIceWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;
    using HoleIceMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct holeiceWeighter : public phys_tools::GenericWeighter<holeiceWeighter<Event,DataType>>{
            private:
                const HoleIceMapType& hole_ice_map_;
                DataType hole_ice_forward_;
                const golemfit::FluxComponent flux_component_;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const HoleIceMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                holeiceWeighter(const HoleIceMapType& hole_ice_map_, DataType hole_ice_forward_, golemfit::FluxComponent flux_component_):
                    hole_ice_map_(hole_ice_map_),
                    hole_ice_forward_(hole_ice_forward_),
                    flux_component_(flux_component_),
                    enforce_needed((not hole_ice_map_.empty()) and IsComponentInMap(hole_ice_map_,flux_component_))
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                    if(hole_ice_map_.empty() or not enforce_needed)
                      return 1.;

                    double cache;
                    if(flux_component_ == golemfit::FluxComponent::atmConv)
                      cache = e.cachedHoleIceConv;
                    else if (flux_component_ == golemfit::FluxComponent::atmPrompt)
                      cache = e.cachedHoleIcePrompt;
                    else if (flux_component_ == golemfit::FluxComponent::diffuseAstro_mu)
                      cache = e.cachedHoleIceAstro;
                    else
                      throw std::runtime_error("HoleIce correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");

                    double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_};

                    auto holeicecorrection = hole_ice_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(holeicecorrection == hole_ice_map_.end())
                        throw std::runtime_error("HoleIce correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    double rate=(*(*holeicecorrection).second)(coordinates);

                    if(cache == 0 and rate > -1.e5)
                      throw std::runtime_error("Badly handled cow.");

                    if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                      // this is a cow, ignore the event
                      return 0.0;
                    }

                    return(pow(10.,rate-cache));
                }
        };

    template<typename Event,unsigned int Dim>
        struct holeiceWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<holeiceWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const HoleIceMapType& hole_ice_map_;
                phys_tools::autodiff::FD<Dim> hole_ice_forward_;
                const golemfit::FluxComponent flux_component_;
                unsigned int didx;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const HoleIceMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                holeiceWeighter(const HoleIceMapType& hole_ice_map_, phys_tools::autodiff::FD<Dim> hole_ice_forward_, golemfit::FluxComponent flux_component_):
                    hole_ice_map_(hole_ice_map_),
                    hole_ice_forward_(hole_ice_forward_),
                    flux_component_(flux_component_),
                    enforce_needed((not hole_ice_map_.empty()) and IsComponentInMap(hole_ice_map_,flux_component_)) {
                    //const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(hole_ice_forward_);
                    for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                        if(hole_ice_forward_.derivative(i)!=0){
                            didx=i;
                            break;
                        }
                    }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                    if(hole_ice_map_.empty() or not enforce_needed)
                      return result_type(1.);

                    double cache;
                    if(flux_component_ == golemfit::FluxComponent::atmConv)
                      cache = e.cachedHoleIceConv;
                    else if (flux_component_ == golemfit::FluxComponent::atmPrompt)
                      cache = e.cachedHoleIcePrompt;
                    else if (flux_component_ == golemfit::FluxComponent::diffuseAstro_mu)
                      cache = e.cachedHoleIceAstro;
                    else
                      throw std::runtime_error("HoleIce correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");

                    double rate, derivative;
                    double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_.value()};
                    int centers[3];

                    auto holeicecorrection = hole_ice_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(holeicecorrection == hole_ice_map_.end())
                        throw std::runtime_error("HoleIce correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    rate=(*(*holeicecorrection).second)(coordinates);

                    if(not (*(*holeicecorrection).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: holeice spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(hole_ice_forward_.value()) + ").");
                    derivative=(*(*holeicecorrection).second).ndsplineeval(coordinates,centers,1<<2);
                    derivative*=hole_ice_forward_.derivative(didx);

                    if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                      // this is a cow, ignore the event
                      result_type r(0.0);
                      r.setDerivative(didx,0.0);
                      return r;
                    }

                    result_type r(rate);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r-cache));
                }
        };

    //used for initializing the per-event dom-eff related caches
    template<typename Event>
        struct HoleIceSetter{
            const HoleIceMapType& hole_ice_map_;
            const double hole_ice_forward_;
            HoleIceSetter(const HoleIceMapType& hole_ice_map_, double hole_ice_forward_):
                hole_ice_map_(hole_ice_map_),hole_ice_forward_(hole_ice_forward_){}

            void setCache(Event& e) const{
                double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_};
                assert(not hole_ice_map_.empty());
                for(auto element : hole_ice_map_){
                  if(element.first.second == static_cast<golemfit::Topology>(e.topology)){
                    double cache = (*element.second)(coordinates);
                    assert(std::isfinite(cache));
                    if(cache == 0){
                      // the cache if the log of the rate; this cases the event to be ignored.
                      cache = -std::numeric_limits<float>::max();
                    }
                    if(element.first.first == golemfit::FluxComponent::atmConv)
                      e.cachedHoleIceConv = cache;
                    else if(element.first.first == golemfit::FluxComponent::atmPrompt)
                      e.cachedHoleIcePrompt = cache;
                    else if(element.first.first == golemfit::FluxComponent::diffuseAstro_mu)
                      e.cachedHoleIceAstro= cache;
                    else
                      throw std::runtime_error("Unimplemented HoleIce component");
                  }
                }
            }
        };

} // close namespace HoleIceWeighterSpace

//================================================================================
// ANISOTROPY WEIGHTER
//================================================================================

namespace AnisotropyWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;
    using AnisotropyMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct anisotropyWeighter : public phys_tools::GenericWeighter<anisotropyWeighter<Event,DataType>>{
            private:
                const AnisotropyMapType& anisotropy_map_;
                DataType anisotropy_scale_;
                const golemfit::FluxComponent flux_component_;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const AnisotropyMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                anisotropyWeighter(const AnisotropyMapType& anisotropy_map_, DataType anisotropy_scale_, golemfit::FluxComponent flux_component_):
                    anisotropy_map_(anisotropy_map_),
                    anisotropy_scale_(anisotropy_scale_),
                    flux_component_(flux_component_),
                    enforce_needed((not anisotropy_map_.empty()) and IsComponentInMap(anisotropy_map_,flux_component_))
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                    if(anisotropy_map_.empty() or not enforce_needed)
                      return 1.;
                    auto cachecorrection = e.cachedAnisotropy.find(flux_component_);
                    if(cachecorrection == e.cachedAnisotropy.end())
                        throw std::runtime_error("Anisotropy correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    float cache = (*cachecorrection).second;
                    double coordinates[2]={log10(e.energy),anisotropy_scale_};
                    if(e.topology == 2) {
                        coordinates[0] = log10(e.length);
                        if(coordinates[0] > 2.3) {
                            return 1.0;
                        }
                    }

                    auto correction = anisotropy_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(correction == anisotropy_map_.end())
                        throw std::runtime_error("Anisotropy correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    double rate=(*(*correction).second)(coordinates);
                    return(pow(10.,rate-cache));
                }
        };

    template<typename Event,unsigned int Dim>
        struct anisotropyWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<anisotropyWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const AnisotropyMapType& anisotropy_map_;
                phys_tools::autodiff::FD<Dim> anisotropy_scale_;
                const golemfit::FluxComponent flux_component_;
                unsigned int didx;
                const bool enforce_needed;
            private:
                bool IsComponentInMap(const AnisotropyMapType& mapa, golemfit::FluxComponent componente){
                  for(auto elo : mapa){
                    if(elo.first.first == componente) return true;
                  }
                  return false;
                }
            public:
                anisotropyWeighter(const AnisotropyMapType& anisotropy_map_, phys_tools::autodiff::FD<Dim> anisotropy_scale_, golemfit::FluxComponent flux_component_):
                    anisotropy_map_(anisotropy_map_),
                    anisotropy_scale_(anisotropy_scale_),
                    flux_component_(flux_component_),
                    enforce_needed((not anisotropy_map_.empty()) and IsComponentInMap(anisotropy_map_,flux_component_)) {
                    //const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(anisotropy_scale_);
                    for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                        if(anisotropy_scale_.derivative(i)!=0){
                            didx=i;
                            break;
                        }
                    }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                    if(anisotropy_map_.empty() or not enforce_needed)
                      return result_type(1.);
                    auto cachecorrection = e.cachedAnisotropy.find(flux_component_);
                    if(cachecorrection == e.cachedAnisotropy.end())
                        throw std::runtime_error("Anisotropy correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    float cache = (*cachecorrection).second;

                    double rate, derivative;
                    double coordinates[2]={log10(e.energy),anisotropy_scale_.value()};
                    if(e.topology == 2) {
                        coordinates[0] = log10(e.length);
                        if(coordinates[0] > 2.3) {
                            return result_type(1.0);
                        }
                    }

                    int centers[2];

                    auto correction = anisotropy_map_.find(std::make_pair(flux_component_,static_cast<golemfit::Topology>(e.topology)));
                    if(correction == anisotropy_map_.end())
                        throw std::runtime_error("Anisotropy correction for " + golemfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + golemfit::GetTopologyName(static_cast<golemfit::Topology>(e.topology)) + " not found.");
                    rate=(*(*correction).second)(coordinates);
                    if(not (*(*correction).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: anisotropy spline topology(" + std::to_string(e.topology) + ") (" + std::to_string(log10(e.length)) + ", " + std::to_string(anisotropy_scale_.value()) + ").");
                    derivative=(*(*correction).second).ndsplineeval(coordinates,centers,1<<1);
                    derivative*=anisotropy_scale_.derivative(didx);

                    result_type r(rate);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r-cache));
                }
        };

    //used for initializing the per-event dom-eff related caches
    template<typename Event>
        struct AnisotropySetter{
            const AnisotropyMapType& anisotropy_map_;
            const double anisotropy_scale_;
            AnisotropySetter(const AnisotropyMapType& anisotropy_map_, double anisotropy_scale_):
                anisotropy_map_(anisotropy_map_),anisotropy_scale_(anisotropy_scale_){}

            void setCache(Event& e) const{
                double coordinates[2]={log10(e.energy),anisotropy_scale_};
                if(e.topology == 2) {
                    coordinates[0] = log10(e.length);
                }
                e.cachedAnisotropy.clear();
                assert(not anisotropy_map_.empty());
                for(auto element : anisotropy_map_){
                  if(element.first.second == static_cast<golemfit::Topology>(e.topology)){
                    if(e.topology == 2 && coordinates[0] > 2.3)
                        e.cachedAnisotropy.insert({element.first.first, 1.0});
                    else
                        e.cachedAnisotropy.insert({element.first.first,(*element.second)(coordinates)});

                  }
                }
                assert(not (e.cachedAnisotropy.empty()));
            }
        };

} // close namespace anisotropyWeighterSpace

//================================================================================
// ATMOSPHERIC DENSITY UNCERATAINTY WEIGHTER
//================================================================================

namespace AtmosphericUncertaintyWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;

    template<typename Event, typename DataType>
        struct atmosphericDensityUncertaintyWeighter : public phys_tools::GenericWeighter<atmosphericDensityUncertaintyWeighter<Event,DataType>>{
            private:
                const std::shared_ptr<splinetable<>> atmospheric_density_uncertainty_spline_;
                DataType scale_;
            public:
                atmosphericDensityUncertaintyWeighter(std::shared_ptr<splinetable<>> atmospheric_density_uncertainty_spline_, DataType scale_):
                    atmospheric_density_uncertainty_spline_(atmospheric_density_uncertainty_spline_),
                    scale_(scale_)
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                  if(atmospheric_density_uncertainty_spline_ == nullptr)
                    return 1.;
                  double coordinates[2]={log10(e.primaryEnergy),cos(e.primaryZenith)};
                  double one_sigma_shift = (*atmospheric_density_uncertainty_spline_)(coordinates);
                  double correction = 1.0+one_sigma_shift*scale_;
                  if(correction<0.0)
                    throw std::runtime_error("GolemFit::AtmosphericDensityUncertaintyWeighter making a neg weight.");
                  return correction;
                }
        };

    template<typename Event,unsigned int Dim>
        struct atmosphericDensityUncertaintyWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<atmosphericDensityUncertaintyWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const std::shared_ptr<splinetable<>> atmospheric_density_uncertainty_spline_;
                phys_tools::autodiff::FD<Dim> scale_;
                unsigned int didx;
            public:
                atmosphericDensityUncertaintyWeighter(std::shared_ptr<splinetable<>> atmospheric_density_uncertainty_spline_, phys_tools::autodiff::FD<Dim> scale_):
                    atmospheric_density_uncertainty_spline_(atmospheric_density_uncertainty_spline_),
                    scale_(scale_) {
                  for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(scale_.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                  }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                  if(atmospheric_density_uncertainty_spline_ == nullptr)
                    return result_type(1.0);
                  double coordinates[2]={log10(e.primaryEnergy),cos(e.primaryZenith)};
                  double one_sigma_shift = (*atmospheric_density_uncertainty_spline_)(coordinates);
                  double rate = 1.0+one_sigma_shift*scale_.value();
                  double derivative = one_sigma_shift;
                  result_type r(rate);
                  r.setDerivative(didx,derivative);
                  return(r);
                }
        };
}// close AtmosphericUncertaintyWeighterSpace

//================================================================================
// KAON ENERGY LOSSES UNCERATAINTY WEIGHTER
//================================================================================

namespace KaonLossesUncertaintyWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;

    template<typename Event, typename DataType>
        struct kaonLossesUncertaintyWeighter: public phys_tools::GenericWeighter<kaonLossesUncertaintyWeighter<Event,DataType>>{
            private:
                const std::shared_ptr<splinetable<>> uncertainty_spline_;
                DataType scale_;
            public:
                kaonLossesUncertaintyWeighter(std::shared_ptr<splinetable<>> uncertainty_spline_, DataType scale_):
                    uncertainty_spline_(uncertainty_spline_),
                    scale_(scale_)
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                  if(uncertainty_spline_ == nullptr){
                    std::cout << "no spline" << std::endl;
                    return 1.;
                  }
                  double coordinates[2]={log10(e.primaryEnergy),cos(e.primaryZenith)};
                  double one_sigma_shift = (*uncertainty_spline_)(coordinates);
                  double correction = 1.0+one_sigma_shift*scale_;
                  if(correction<0.0)
                    throw std::runtime_error("GolemFit::AtmosphericDensityUncertaintyWeighter making a neg weight.");
                  return correction;
                }
        };

    template<typename Event,unsigned int Dim>
        struct kaonLossesUncertaintyWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<kaonLossesUncertaintyWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const std::shared_ptr<splinetable<>> uncertainty_spline_;
                phys_tools::autodiff::FD<Dim> scale_;
                unsigned int didx;
            public:
                kaonLossesUncertaintyWeighter(std::shared_ptr<splinetable<>> uncertainty_spline_, phys_tools::autodiff::FD<Dim> scale_):
                    uncertainty_spline_(uncertainty_spline_),
                    scale_(scale_) {
                  for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(scale_.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                  }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                  if(uncertainty_spline_ == nullptr)
                    return result_type(1.0);
                  double coordinates[2]={log10(e.primaryEnergy),cos(e.primaryZenith)};
                  double one_sigma_shift = (*uncertainty_spline_)(coordinates);
                  double rate = 1.0+one_sigma_shift*scale_.value();
                  double derivative = one_sigma_shift;
                  result_type r(rate);
                  r.setDerivative(didx,derivative);
                  return(r);
                }
        };
}// close KaonLossesUncertaintyWeighterSpace


//================================================================================
// ICE GRADIENT UNCERATAINTY WEIGHTER
//================================================================================

namespace IceGradientsWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;

    template<typename Event, typename DataType>
        struct icegradientWeighter: public phys_tools::GenericWeighter<icegradientWeighter<Event,DataType>>{
            private:
                const std::vector<double> bin_edges_;
                //const nusquids::marray<double,1> gradient_;
                const std::vector<double> gradient_;
                DataType scale_;
                bool gradients_loaded_;
            public:
                icegradientWeighter():gradients_loaded_(false){}

                //icegradientWeighter(std::vector<double> bin_edges_, nusquids::marray<double,1> gradient_, DataType scale_):
                icegradientWeighter(std::vector<double> bin_edges_, std::vector<double> gradient_, DataType scale_):
                    bin_edges_(bin_edges_),
                    gradient_(gradient_),
                    scale_(scale_),
                    gradients_loaded_(true)
                {}

                size_t Get_index(double energy) const {
                  //for(auto edge: bin_edges_)
                  //  std::cout << edge << " ";
                  //std::cout << energy << std::endl;
                  auto bedgeit=std::lower_bound(bin_edges_.begin(),bin_edges_.end(),log10(energy));
                  if(bedgeit==bin_edges_.end())
                    return std::numeric_limits<size_t>::max();
                    //throw std::runtime_error("CA:End of the array");
                  if(bedgeit!=bin_edges_.begin())
                    bedgeit--;
                  size_t index_M=std::distance(bin_edges_.begin(),bedgeit);
                  return index_M;
                }

                using result_type=double;
                result_type operator()(const Event& e) const{
                  if(not gradients_loaded_) return 1.;
                  if(gradient_.size() == 0)
                    return 1.;
                  size_t index = Get_index(e.energy);
                  //std::cout << index << " " <<  gradient_.size() << " " << scale_ << std::endl;
                  if(index>=bin_edges_.size())
                    return 1.;
                  double rel = (1. + gradient_[index]*scale_);
                  if(rel<0)
                    throw std::runtime_error("CA:Negative weight");
                  return rel;
                }

                icegradientWeighter(const icegradientWeighter& other):
                  bin_edges_((other.bin_edges_)),
                  gradient_((other.gradient_)),
                  scale_((other.scale_)),
                  gradients_loaded_((other.gradients_loaded_))
                {}

                icegradientWeighter(icegradientWeighter&& other):
                  bin_edges_(std::move(other.bin_edges_)),
                  gradient_(std::move(other.gradient_)),
                  scale_(std::move(other.scale_)),
                  gradients_loaded_(std::move(other.gradients_loaded_))
                {
                  other.gradients_loaded_ = false;
                }

                icegradientWeighter& operator=(const icegradientWeighter& right)
                {
                  // some evilness.
                  if (this == &right) return *this;
                  this->~icegradientWeighter();
                  new (this) icegradientWeighter(right);
                  return *this;
                }
        };

    template<typename Event,unsigned int Dim>
        struct icegradientWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<icegradientWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const std::vector<double> bin_edges_;
                //`const nusquids::marray<double,1> gradient_;
                const std::vector<double> gradient_;
                phys_tools::autodiff::FD<Dim> scale_;
                bool gradients_loaded_;
                unsigned int didx;
            public:
                icegradientWeighter():gradients_loaded_(false){}
                //icegradientWeighter(std::vector<double> bin_edges_, nusquids::marray<double,1> gradient_, phys_tools::autodiff::FD<Dim> scale_):
                icegradientWeighter(std::vector<double> bin_edges_, std::vector<double> gradient_, phys_tools::autodiff::FD<Dim> scale_):
                    bin_edges_(bin_edges_),
                    gradient_(gradient_),
                    scale_(scale_),
                    gradients_loaded_(true)
                {
                    for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                      if(scale_.derivative(i)!=0){
                          didx=i;
                          break;
                      }
                    }
                }

                size_t Get_index(double energy) const {
                  auto bedgeit=std::lower_bound(bin_edges_.begin(),bin_edges_.end(),log10(energy));
                  if(bedgeit==bin_edges_.end())
                    return std::numeric_limits<size_t>::max();
                    //throw std::runtime_error("CA:End of the array");
                  if(bedgeit!=bin_edges_.begin())
                    bedgeit--;
                  size_t index_M=std::distance(bin_edges_.begin(),bedgeit);
                  return index_M;
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const {
                  if(not gradients_loaded_) return result_type(1.);
                  if(gradient_.size() == 0)
                    return result_type(1.0);
                  size_t index = Get_index(e.energy);
                  if(index>=bin_edges_.size())
                    return result_type(1.0);
                  double derivative = gradient_[index];
                  double rel = 1+derivative*scale_.value();
                  result_type r(rel);
                  r.setDerivative(didx,derivative);
                  return(r);
                }

                icegradientWeighter(const icegradientWeighter& other):
                  bin_edges_((other.bin_edges_)),
                  gradient_((other.gradient_)),
                  scale_((other.scale_)),
                  gradients_loaded_((other.gradients_loaded_)),
                  didx((other.didx))
                {}

                icegradientWeighter(icegradientWeighter&& other):
                  bin_edges_(std::move(other.bin_edges_)),
                  gradient_(std::move(other.gradient_)),
                  scale_(std::move(other.scale_)),
                  gradients_loaded_(std::move(other.gradients_loaded_)),
                  didx(std::move(other.didx))
                {
                  other.gradients_loaded_ = false;
                }

                icegradientWeighter& operator=(const icegradientWeighter& right)
                {
                  // some evilness.
                  if (this == &right) return *this;
                  this->~icegradientWeighter();
                  new (this) icegradientWeighter(right);
                  return *this;
                }
        };

}// close IceGradientsWeighterSpace

//================================================================================
// ASTRO MODEL WEIGHTER
//================================================================================

namespace AstroModelWeighterSpaceLegacy {
    using namespace photospline;
    using namespace golemfit;
    using namespace phys_tools::autodiff;
    using AstroMapType=std::map<AstrophysicalNeutrinoModel,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct astroModelWeighter : public phys_tools::GenericWeighter<astroModelWeighter<Event,DataType>>{
            private:
                const AstroMapType& astroMap;
                AstrophysicalNeutrinoModel model;
                DataType loge_shift;
            public:
                astroModelWeighter(AstrophysicalNeutrinoModel model, const AstroMapType& astroMap, DataType loge_shift):
                    astroMap(astroMap),model(model),loge_shift(loge_shift){}

                using result_type=double;
                result_type operator()(const Event& e) const{
                    double coordinates[1]={log10(e.primaryEnergy) - loge_shift};

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    double log10flux =(*(*spline).second)(coordinates);
                    return(pow(10.,log10flux));
                }
        };

    template<typename Event,unsigned int Dim>
        struct astroModelWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<astroModelWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                AstroMapType astroMap;
                AstrophysicalNeutrinoModel model;
                phys_tools::autodiff::FD<Dim> loge_shift;
                unsigned int didx;
            public:
                astroModelWeighter(AstrophysicalNeutrinoModel model, AstroMapType astroMap, phys_tools::autodiff::FD<Dim> loge_shift):
                    astroMap(astroMap),model(model),loge_shift(loge_shift)
            {
                const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(loge_shift);
                for(unsigned int i=0; i<n; i++){
                    if(loge_shift.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                }
            }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                    double flux, derivative;
                    double coordinates[1]={log10(e.primaryEnergy) - loge_shift.value()};
                    int centers[1];

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    flux =(*(*spline).second)(coordinates);
                    if(not (*(*spline).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: logE " + std::to_string(coordinates[0]) + " with shift " + std::to_string(loge_shift.value()));
                    derivative=(-1.)*(*(*spline).second).ndsplineeval(coordinates,centers,1<<0);
                    derivative*=loge_shift.derivative(didx);

                    result_type r(flux);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r));
                }
        };
} // close namespace astroweighter legacy namespace

//================================================================================
// ASTRO MODEL WEIGHTER
//================================================================================

namespace AstroModelWeighterSpace {
    using namespace photospline;
    using namespace golemfit;
    using namespace phys_tools::autodiff;
    using AstroMapType=std::map<AstrophysicalNeutrinoModel,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct astroModelWeighter : public phys_tools::GenericWeighter<astroModelWeighter<Event,DataType>>{
            private:
                const AstroMapType& astroMap;
                AstrophysicalNeutrinoModel model;
                DataType loge_shift;
                const unsigned int spline_dimension;
                const double minE,maxE;
                double minCosth, maxCosth;
            public:
                astroModelWeighter(AstrophysicalNeutrinoModel model, const AstroMapType& astroMap, DataType loge_shift):
                    astroMap(astroMap),model(model),loge_shift(loge_shift),
                    spline_dimension((*astroMap.find(model)).second->get_ndim()),
                    minE((*astroMap.find(model)).second->lower_extent(0)),
                    maxE((*astroMap.find(model)).second->upper_extent(0)){
                      if(spline_dimension > 1){
                        minCosth = (*astroMap.find(model)).second->lower_extent(1);
                        maxCosth = (*astroMap.find(model)).second->upper_extent(1);
                      }
                    }

                using result_type=double;
                result_type operator()(const Event& e) const{
                  if(spline_dimension == 1){
                    double coordinates[1]={log10(e.primaryEnergy) - loge_shift};
                    if(coordinates[0] >= maxE or coordinates[0] <= minE)
                      return 0.0;

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    double log10flux =(*(*spline).second)(coordinates);
                    return(pow(10.,log10flux));
                  } else if (spline_dimension == 2) {
                    double coordinates[2]={log10(e.primaryEnergy) - loge_shift, cos(e.primaryZenith)};
                    if(coordinates[0] >= maxE or coordinates[0] <= minE)
                      return 0.0;
                    if(coordinates[1] >= maxCosth or coordinates[1] <= minCosth)
                      return 0.0;

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    double log10flux =(*(*spline).second)(coordinates);
                    return(pow(10.,log10flux));
                  } else
                    throw std::runtime_error("astroModelWeighter::operator(): invalid spline dimension");
                }
        };

    template<typename Event,unsigned int Dim>
        struct astroModelWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<astroModelWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                AstroMapType astroMap;
                AstrophysicalNeutrinoModel model;
                phys_tools::autodiff::FD<Dim> loge_shift;
                unsigned int didx;
                const unsigned int spline_dimension;
                const double minE,maxE;
                double minCosth, maxCosth;
            public:
                astroModelWeighter(AstrophysicalNeutrinoModel model, AstroMapType astroMap, phys_tools::autodiff::FD<Dim> loge_shift):
                    astroMap(astroMap),model(model),loge_shift(loge_shift),
                    spline_dimension((*astroMap.find(model)).second->get_ndim()),
                    minE((*astroMap.find(model)).second->lower_extent(0)),
                    maxE((*astroMap.find(model)).second->upper_extent(0))
            {
                const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(loge_shift);
                for(unsigned int i=0; i<n; i++){
                    if(loge_shift.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                }

                if(spline_dimension > 1){
                  minCosth = (*astroMap.find(model)).second->lower_extent(1);
                  maxCosth = (*astroMap.find(model)).second->upper_extent(1);
                }
            }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                  if(spline_dimension == 1){
                    double flux, derivative;
                    double coordinates[1]={log10(e.primaryEnergy) - loge_shift.value()};
                    int centers[1];

                    if(coordinates[0] >= maxE or coordinates[0] <= minE)
                      return result_type(0.0);

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    flux =(*(*spline).second)(coordinates);
                    if(not (*(*spline).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: logE " + std::to_string(coordinates[0]) + " with shift " + std::to_string(loge_shift.value()));
                    derivative=(-1.)*(*(*spline).second).ndsplineeval(coordinates,centers,1<<0);
                    derivative*=loge_shift.derivative(didx);

                    result_type r(flux);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r));
                  } else if (spline_dimension == 2) {
                    // todo check derivative with numerical calculation
                    double flux, derivative;
                    double coordinates[2]={log10(e.primaryEnergy) - loge_shift.value(), cos(e.primaryZenith)};
                    int centers[2];

                    if(coordinates[0] >= maxE or coordinates[0] <= minE)
                      return result_type(0.0);
                    if(coordinates[1] >= maxCosth or coordinates[1] <= minCosth)
                      return result_type(0.0);

                    auto spline = astroMap.find(model);
                    if( spline == astroMap.end() )
                        throw std::runtime_error("Astrophysical model not found.");
                    flux =(*(*spline).second)(coordinates);
                    if(not (*(*spline).second).searchcenters(coordinates,centers))
                        throw std::runtime_error("Out of phase space: logE " + std::to_string(coordinates[0]) + " with shift " + std::to_string(loge_shift.value()) + " and costh " + std::to_string(coordinates[1]));
                    derivative=(-1.)*(*(*spline).second).ndsplineeval(coordinates,centers,1<<0);
                    derivative*=loge_shift.derivative(didx);

                    result_type r(flux);
                    r.setDerivative(didx,derivative);
                    return(pow(10.,r));
                  } else
                    throw std::runtime_error("astroModelWeighter::operator(): invalid spline dimension");
                }
        };
} // close namespace astroweighter namespace


//================================================================================
// Attenuation WEIGHTER
//================================================================================

namespace DirectAttenuationWeighterSpace {

    using namespace nufate;
    using namespace nusquids;
    using namespace photospline;
    using namespace golemfit;
    using namespace phys_tools::autodiff;
    using AstroMapType=std::map<AstrophysicalNeutrinoModel,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct attenuationWeighter : public phys_tools::GenericWeighter<attenuationWeighter<Event,DataType>>{
            private:
                FluxComponent flux_component;
                DataType xs_scaling;
                //DataType nu_xs_scaling;
                //DataType nubar_xs_scaling;
                std::map<LW::ParticleType,nuFATE> nufate_baseline, nufate;
                const std::string xs_path;
                const bool secondaries_included;
                std::vector<double> costh_range;
                std::map<LW::ParticleType,std::vector<AkimaSpline>> get_attenuation_change;

                double SpitAproxSpectrum(FluxComponent flux_component) const {
                   switch(flux_component){
                     case FluxComponent::atmConv : return -3.75; break;
                     case FluxComponent::atmPrompt : return -2.75; break;
                     case FluxComponent::diffuseAstro_e  : return -2.5; break;
                     case FluxComponent::diffuseAstro_mu : return -2.5; break;
                     case FluxComponent::diffuseAstro_tau: return -2.5; break;
                     default: assert(false && "Component de flujo no especificada.");
                   }
                }
                int TonuFATE(LW::ParticleType particle_type) const {
                  switch(particle_type){
                    case LW::ParticleType::NuMu: return 2; break;
                    case LW::ParticleType::NuMuBar: return -2; break;
                    default: assert(false && "Traductor de nuFATE no valido.");
                  }
                }
            public:
                attenuationWeighter(FluxComponent flux_component,
                                    DataType xs_scaling,
                                    //DataType nu_xs_scaling, DataType nubar_xs_scaling,
                                    std::string xs_path, bool secondaries_included):
                    flux_component(flux_component),
                    xs_scaling(xs_scaling),
                    //nu_xs_scaling(nu_xs_scaling), nubar_xs_scaling(nubar_xs_scaling),
                    xs_path(xs_path), secondaries_included(secondaries_included) {
                      double gamma = SpitAproxSpectrum(flux_component);
                      nufate_baseline = {{LW::ParticleType::NuMu,nuFATE(2,gamma, xs_path, secondaries_included)},
                                         {LW::ParticleType::NuMuBar,nuFATE(-2,gamma, xs_path, secondaries_included)}};

                      for(auto primary_type : { LW::ParticleType::NuMu, LW::ParticleType::NuMuBar }) {
                        int ptype = TonuFATE(primary_type);

                        std::vector<double> energy_nodes = nufate_baseline.at(primary_type).getEnergyNodes();
                        std::vector<double> total_xs = nufate_baseline.at(primary_type).getTotalCrossSections();
                        auto diff_xs_nc = nufate_baseline.at(primary_type).getNCDifferentialCrossSections();

                        double scale_factor = xs_scaling;
                        //double scale_factor = isNeutrino(primary_type) ? nu_xs_scaling : nubar_xs_scaling;

                        std::vector<std::vector<double>> dxsnc;
                        unsigned int NumNodes = diff_xs_nc.n_;
                        std::shared_ptr<double> XSVec = diff_xs_nc.vec_;
                        for(unsigned int i = 0; i<NumNodes; i++){
                          std::vector<double> dsigma_row(NumNodes);
                          for(unsigned int j=0; j<NumNodes; j++){
                            dsigma_row[j] = *(XSVec.get()+i*NumNodes+j)*scale_factor;
                          }
                          dxsnc.push_back(dsigma_row);
                          total_xs[i] = total_xs[i]*scale_factor;
                        }

                        nufate.insert({primary_type,nuFATE(ptype,gamma,energy_nodes,total_xs,dxsnc,secondaries_included)});
                      }

                      double costh_min = -1.;
                      double costh_max =  1.;
                      double costh_delta = 0.01;
                      costh_range.clear();
                      for(double costh = costh_min; costh <= costh_max; costh += costh_delta){
                        costh_range.push_back(costh);
                      }

                      for(auto primary_type : { LW::ParticleType::NuMu, LW::ParticleType::NuMuBar }) {
                        auto baseline_nufate = nufate_baseline.at(primary_type);
                        auto alternative_nufate = nufate.at(primary_type);

                        std::vector<double> energy_nodes = baseline_nufate.getEnergyNodes();

                        std::vector<AkimaSpline> akima_slice;
                        for(double costh : costh_range){
                          double number_of_targets = baseline_nufate.getEarthColumnDensity(acos(costh));
                          std::vector<double> baseline_attenuation = baseline_nufate.getRelativeAttenuation(number_of_targets);
                          std::vector<double> alternative_attenuation = alternative_nufate.getRelativeAttenuation(number_of_targets);
                          std::vector<double> ratio_of_attenuations;
                          for(unsigned int i = 0; i < baseline_attenuation.size(); i++){
                            ratio_of_attenuations.push_back(alternative_attenuation[i]/baseline_attenuation[i]);
                          }
                          akima_slice.push_back(AkimaSpline(energy_nodes, ratio_of_attenuations));
                        }
                        get_attenuation_change.insert({primary_type, akima_slice});
                      }
                }

                size_t Get_index(double costh) const {
                  //for(auto edge: bin_edges_)
                  //  std::cout << edge << " ";
                  //std::cout << costh << std::endl;
                  auto bedgeit=std::lower_bound(costh_range.begin(),costh_range.end(),costh);
                  if(bedgeit==costh_range.end())
                    return std::numeric_limits<size_t>::max();
                    //throw std::runtime_error("CA:End of the array");
                  if(bedgeit!=costh_range.begin())
                    bedgeit--;
                  size_t index_M=std::distance(costh_range.begin(),bedgeit);
                  return index_M;
                }

                using result_type=double;
                result_type operator()(const Event& e) const{
                  double costh = cos(e.primaryZenith);
                  return get_attenuation_change.at(e.primaryType)[Get_index(costh)](e.primaryEnergy);
                }
        };

    template<typename Event,unsigned int Dim>
        struct attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                FluxComponent flux_component;
                std::map<LW::ParticleType,nuFATE> nufate_baseline, nufate;
                const std::string xs_path;
                const bool secondaries_included;
                std::vector<double> costh_range;
                std::map<LW::ParticleType,std::vector<AkimaSpline>> get_attenuation_change;
                phys_tools::autodiff::FD<Dim> xs_scaling;
                //DataType nu_xs_scaling;
                //DataType nubar_xs_scaling;
                unsigned int didx;
            public:
                attenuationWeighter(FluxComponent flux_component,
                                    phys_tools::autodiff::FD<Dim> xs_scaling,
                                    std::string xs_path, bool secondaries_included):
                    flux_component(flux_component),
                    xs_scaling(xs_scaling),
                    //nu_xs_scaling(nu_xs_scaling), nubar_xs_scaling(nubar_xs_scaling),
                    xs_path(xs_path), secondaries_included(secondaries_included)
                {
                  const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(xs_scaling);
                  for(unsigned int i=0; i<n; i++){
                      if(xs_scaling.derivative(i)!=0){
                          didx=i;
                          break;
                      }
                  }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                  return result_type(0.0);
                }
        };

} // close namespace earth absorption direct reweighter

namespace SplineAttenuationWeighterSpace {
    using namespace photospline;
    using namespace phys_tools::autodiff;
    using AttenuationMapType=std::map<std::pair<golemfit::FluxComponent,LW::ParticleType>,std::shared_ptr<splinetable<>>>;

    template<typename Event, typename DataType>
        struct attenuationWeighter : public phys_tools::GenericWeighter<attenuationWeighter<Event,DataType>>{
            private:
                const AttenuationMapType & attenuation_uncertainty_spline_;
                golemfit::FluxComponent flux_component_;
                DataType scale_nu_;
                DataType scale_nubar_;
            public:
                attenuationWeighter(const AttenuationMapType &  attenuation_uncertainty_spline_, golemfit::FluxComponent flux_component_, DataType scale_nu_, DataType scale_nubar_):
                    attenuation_uncertainty_spline_(attenuation_uncertainty_spline_),
                    flux_component_(flux_component_),
                    scale_nu_(scale_nu_),
                    scale_nubar_(scale_nubar_)
                {}

                using result_type=double;
                result_type operator()(const Event& e) const{
                  auto attenuation_ratio_spline = attenuation_uncertainty_spline_.find(std::make_pair(flux_component_,e.primaryType));
                  if(attenuation_ratio_spline == attenuation_uncertainty_spline_.end())
                     return result_type(1.0);
                  DataType scale;
                  if(isNeutrino(e.primaryType))
                    scale = scale_nu_;
                  else
                    scale = scale_nubar_;
                  double coordinates[3]={log10(e.primaryEnergy),cos(e.primaryZenith), scale};
                  double correction = (*(*attenuation_ratio_spline).second)(coordinates);
                  if(correction<0.0)
                    throw std::runtime_error("GolemFit::attenuation weighter making a neg weight.");
                  return correction;
                }
        };

    template<typename Event,unsigned int Dim>
        struct attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
            private:
                const AttenuationMapType &  attenuation_uncertainty_spline_;
                golemfit::FluxComponent flux_component_;
                phys_tools::autodiff::FD<Dim> scale_nu_;
                phys_tools::autodiff::FD<Dim> scale_nubar_;
                unsigned int didx_nu_,didx_nubar_;
            public:
                attenuationWeighter(const AttenuationMapType &  attenuation_uncertainty_spline_, golemfit::FluxComponent flux_component_,
                    phys_tools::autodiff::FD<Dim> scale_nu_,
                    phys_tools::autodiff::FD<Dim> scale_nubar_):
                    attenuation_uncertainty_spline_(attenuation_uncertainty_spline_),
                    flux_component_(flux_component_),
                    scale_nu_(scale_nu_),
                    scale_nubar_(scale_nubar_)
                {
                  for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(scale_nu_.derivative(i)!=0){
                        didx_nu_=i;
                        break;
                    }
                  }
                  for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(scale_nubar_.derivative(i)!=0){
                        didx_nubar_=i;
                        break;
                    }
                  }
                }

                using result_type=phys_tools::autodiff::FD<Dim>;
                result_type operator()(const Event& e) const{
                  auto attenuation_ratio_spline = attenuation_uncertainty_spline_.find(std::make_pair(flux_component_,e.primaryType));
                  if(attenuation_ratio_spline == attenuation_uncertainty_spline_.end())
                     return result_type(1.0);
                  double scale;
                  unsigned int didx;
                  if(isNeutrino(e.primaryType)){
                    scale = scale_nu_.value();
                    didx = didx_nu_;
                  } else {
                    scale = scale_nubar_.value();
                    didx = didx_nubar_;
                  }
                  double coordinates[3]={log10(e.primaryEnergy),cos(e.primaryZenith), scale};
                  double correction = (*(*attenuation_ratio_spline).second)(coordinates);

                  int centers[3];
                  if(not (*(*attenuation_ratio_spline).second).searchcenters(coordinates,centers))
                      throw std::runtime_error("Out of phase space: attenuation spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(scale) + ").");
                  double derivative=(*(*attenuation_ratio_spline).second).ndsplineeval(coordinates,centers,1<<2);

                  result_type c(correction);
                  c.setDerivative(didx,derivative);
                  return(c);
                }
        };

} // close namespace earth absorption spline reweighter

using namespace DOMEff5;
using namespace HoleIceWeighterSpace;
using namespace AnisotropyWeighterSpace;
using namespace AstroModelWeighterSpace;
using namespace AtmosphericUncertaintyWeighterSpace;
using namespace KaonLossesUncertaintyWeighterSpace;
using namespace IceGradientsWeighterSpace;
using namespace SplineAttenuationWeighterSpace;

//================================================================================
// THE WEIGHT MAKER
//================================================================================

// This function construct a Weighter object to weight each event independently.
//
// "Time can't be measured in days the way money is measured in pesos and centavos,
// because all pesos are equal, while every day, perhaps every hour, is different."
// JLB

namespace sterile {

struct WeighterMaker{
    private:
        static constexpr double medianConvEnergy=2020;
        //static constexpr double medianPromptEnergy=7887;
        static constexpr double medianPromptEnergy=2020;
        static constexpr double astroPivotEnergy=1.0e5;
        // DOM efficienfy splines;
        using DOMMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;
        using HQDOMMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;
        using HoleIceMapType=std::map<std::pair<golemfit::FluxComponent,golemfit::Topology>,std::shared_ptr<splinetable<>>>;
        using AttenuationMapType=std::map<std::pair<golemfit::FluxComponent, LW::ParticleType>,std::shared_ptr<splinetable<>>>;
        //using IceGradientMapType=std::map<golemfit::IceParameter,std::pair<std::vector<double>,std::vector<double>>;
        using IceGradientMapType=std::pair<std::vector<double>, std::vector<std::vector<double>>>;
        std::shared_ptr<const DOMMapType> domefficiencySplines_;
        std::shared_ptr<const HQDOMMapType> hqdomefficiencySplines_;
        std::shared_ptr<const HoleIceMapType> holeiceSplines_;
        std::shared_ptr<const AttenuationMapType> attenuationSplines_;
        std::shared_ptr<const AstroMapType> astrophysicalNeutrinoModels_;
        std::shared_ptr<splinetable<>> atmosphericDensityUncertaintySpline_;
        std::shared_ptr<splinetable<>> kaonLossesUncertaintySpline_;
        //std::pair<std::vector<double>,std::vector<nusquids::marray<double,2>>> ice_gradients_;
        IceGradientMapType ice_gradients_;
        std::string nufate_xs_path_ = "";
        bool domeff_splines_loaded_ = false;
        bool hqdomeff_splines_loaded_ = false;
        bool astrophysical_models_splines_loaded_ = false;
        bool holeice_splines_loaded_ = false;
        bool attenuation_splines_loaded_ = false;
        bool atmospheric_density_uncertainty_spline_loaded_ = false;
        bool kaon_losses_spline_loaded_ = false;
        bool ice_gradients_loaded_ = false;
        golemfit::DiffuseFitType fit_type;
        golemfit::SteeringParams const * steering;
    public:
        // default constructor // bad bad
        WeighterMaker():
          domefficiencySplines_(nullptr),
          hqdomefficiencySplines_(nullptr),
          holeiceSplines_(nullptr),
          astrophysicalNeutrinoModels_(nullptr),
          atmosphericDensityUncertaintySpline_(nullptr),
          nufate_xs_path_(""),
          domeff_splines_loaded_(false),
          hqdomeff_splines_loaded_(false),
          astrophysical_models_splines_loaded_(false),
          holeice_splines_loaded_(false),
          attenuation_splines_loaded_(false),
          atmospheric_density_uncertainty_spline_loaded_(false),
          kaon_losses_spline_loaded_(false),
          ice_gradients_loaded_(false),
          fit_type(golemfit::DiffuseFitType::SinglePowerLaw),
          steering(nullptr)
        {}
        WeighterMaker(const DOMMapType& dom_splines, const DOMMapType& hqdom_splines,const HoleIceMapType& holeice_splines, const AttenuationMapType& attenuation_splines, const AstroMapType& astro_splines, std::shared_ptr<splinetable<>> atmosphericDensityUncertaintySpline, std::shared_ptr<splinetable<>> kaonLossesUncertaintySpline, const IceGradientMapType ice_gradients, const golemfit::SteeringParams& steeringParams, const golemfit::DataPaths& datapaths):
            domefficiencySplines_(std::make_shared<const DOMMapType>(dom_splines)),
            hqdomefficiencySplines_(std::make_shared<const HQDOMMapType>(hqdom_splines)),
            holeiceSplines_(std::make_shared<const HoleIceMapType>(holeice_splines)),
            attenuationSplines_(std::make_shared<const AttenuationMapType>(attenuation_splines)),
            astrophysicalNeutrinoModels_(std::make_shared<const AstroMapType>(astro_splines)),
            atmosphericDensityUncertaintySpline_(atmosphericDensityUncertaintySpline),
            kaonLossesUncertaintySpline_(kaonLossesUncertaintySpline),
            ice_gradients_(ice_gradients),
            nufate_xs_path_(datapaths.nufate_xs_table_location),
            domeff_splines_loaded_(true),
            hqdomeff_splines_loaded_(true),
            astrophysical_models_splines_loaded_(true),
            holeice_splines_loaded_(true),
            attenuation_splines_loaded_(true),
            atmospheric_density_uncertainty_spline_loaded_(true),
            kaon_losses_spline_loaded_(true),
            ice_gradients_loaded_(true),
            fit_type(steeringParams.diffuse_fit_type),
            steering(&steeringParams)
        {
          if(ice_gradients_.first.size() == 0 or ice_gradients_.second.size() ==0) ice_gradients_loaded_ = false;
        }
        // copy constructor // CA
        WeighterMaker(const WeighterMaker& guy):
          domefficiencySplines_(guy.domefficiencySplines_),
          hqdomefficiencySplines_(guy.hqdomefficiencySplines_),
          holeiceSplines_(guy.holeiceSplines_),
          attenuationSplines_(guy.attenuationSplines_),
          astrophysicalNeutrinoModels_(guy.astrophysicalNeutrinoModels_),
          atmosphericDensityUncertaintySpline_(guy.atmosphericDensityUncertaintySpline_),
          kaonLossesUncertaintySpline_(guy.kaonLossesUncertaintySpline_),
          ice_gradients_(guy.ice_gradients_),
          nufate_xs_path_(guy.nufate_xs_path_),
          domeff_splines_loaded_(guy.domeff_splines_loaded_),
          hqdomeff_splines_loaded_(guy.hqdomeff_splines_loaded_),
          astrophysical_models_splines_loaded_(guy.astrophysical_models_splines_loaded_),
          holeice_splines_loaded_(guy.holeice_splines_loaded_),
          attenuation_splines_loaded_(guy.attenuation_splines_loaded_),
          atmospheric_density_uncertainty_spline_loaded_(guy.atmospheric_density_uncertainty_spline_loaded_),
          kaon_losses_spline_loaded_(guy.kaon_losses_spline_loaded_),
          ice_gradients_loaded_(guy.ice_gradients_loaded_),
          fit_type(guy.fit_type),
          steering(guy.steering)
        {}
        // assign guy // CA
        WeighterMaker& operator=(WeighterMaker& other){
          if(&other==this)
            return(*this);

          domefficiencySplines_ = other.domefficiencySplines_;
          hqdomefficiencySplines_ = other.hqdomefficiencySplines_;
          holeiceSplines_ = other.holeiceSplines_;
          attenuationSplines_ = other.attenuationSplines_;
          astrophysicalNeutrinoModels_ = other.astrophysicalNeutrinoModels_;
          ice_gradients_ = other.ice_gradients_;
          ice_gradients_loaded_ = other.ice_gradients_loaded_;
          nufate_xs_path_ = other.nufate_xs_path_;
          domeff_splines_loaded_ = other.domeff_splines_loaded_;
          hqdomeff_splines_loaded_ = other.hqdomeff_splines_loaded_;
          atmosphericDensityUncertaintySpline_ = other.atmosphericDensityUncertaintySpline_;
          kaonLossesUncertaintySpline_ = other.kaonLossesUncertaintySpline_;
          astrophysical_models_splines_loaded_ = other.astrophysical_models_splines_loaded_;
          holeice_splines_loaded_ = other.holeice_splines_loaded_;
          attenuation_splines_loaded_ = other.attenuation_splines_loaded_;
          atmospheric_density_uncertainty_spline_loaded_ = other.atmospheric_density_uncertainty_spline_loaded_;
          kaon_losses_spline_loaded_ = other.kaon_losses_spline_loaded_;
          fit_type = other.fit_type;
          steering = other.steering;

          return(*this);
        }
        // move guy // CA
        WeighterMaker& operator=(WeighterMaker&& other){
          if(&other==this)
            return(*this);

          domefficiencySplines_ = std::move(other.domefficiencySplines_);
          hqdomefficiencySplines_ = std::move(other.hqdomefficiencySplines_);
          holeiceSplines_ = std::move(other.holeiceSplines_);
          attenuationSplines_ = std::move(other.attenuationSplines_);
          astrophysicalNeutrinoModels_ = std::move(other.astrophysicalNeutrinoModels_);
          atmosphericDensityUncertaintySpline_ = std::move(other.atmosphericDensityUncertaintySpline_);
          kaonLossesUncertaintySpline_ = std::move(other.kaonLossesUncertaintySpline_);
          ice_gradients_ = std::move(other.ice_gradients_);
          ice_gradients_loaded_ = std::move(other.ice_gradients_loaded_);
          domeff_splines_loaded_ = std::move(other.domeff_splines_loaded_);
          nufate_xs_path_ = std::move(other.nufate_xs_path_);
          hqdomeff_splines_loaded_ = std::move(other.hqdomeff_splines_loaded_);
          astrophysical_models_splines_loaded_ = std::move(other.astrophysical_models_splines_loaded_);
          holeice_splines_loaded_ = std::move(other.holeice_splines_loaded_);
          attenuation_splines_loaded_ = std::move(other.attenuation_splines_loaded_);
          atmospheric_density_uncertainty_spline_loaded_ = std::move(other.atmospheric_density_uncertainty_spline_loaded_);
          kaon_losses_spline_loaded_ = std::move(other.kaon_losses_spline_loaded_);
          fit_type = std::move(other.fit_type);
          steering = std::move(other.steering);

          return(*this);
        }
    public:
        void SetDOMEfficiencySplines(const DOMMapType& dom_splines){ domefficiencySplines_ = std::make_shared<const DOMMapType>(dom_splines); domeff_splines_loaded_ = true; }
        void SetHQDOMEfficiencySplines(const HQDOMMapType& hqdom_splines){ hqdomefficiencySplines_ = std::make_shared<const HQDOMMapType>(hqdom_splines); hqdomeff_splines_loaded_ = true; }
        void SetHoleIceSplines(const HoleIceMapType& holeice_splines){ holeiceSplines_ = std::make_shared<const HoleIceMapType>(holeice_splines); holeice_splines_loaded_ = true; }
        void SetAttenuationSplines(const AttenuationMapType& attenuation_splines){ attenuationSplines_ = std::make_shared<const AttenuationMapType>(attenuation_splines); attenuation_splines_loaded_ = true; }
        void SetAstrophysicalModelSplines(const AstroMapType& astro_splines){ astrophysicalNeutrinoModels_ = std::make_shared<const AstroMapType>(astro_splines); astrophysical_models_splines_loaded_ = true; }
        void SetSteering(const SteeringParams& steeringParams){steering = &steeringParams;}
        void SetAtmosphericDensityUncertaintySpline(std::shared_ptr<splinetable<>> atmosphericDensityUncertaintySpline){
          atmosphericDensityUncertaintySpline_ = atmosphericDensityUncertaintySpline;
          atmospheric_density_uncertainty_spline_loaded_ = true; }
        void SetKaonEnergyLossesUncertaintySpline(std::shared_ptr<splinetable<>> kaonLossesUncertaintySpline){
          kaonLossesUncertaintySpline_ = kaonLossesUncertaintySpline;
          kaon_losses_spline_loaded_= true; }
        void SetIceGradients(std::pair<std::vector<double>,std::vector<std::vector<double>>> ice_Gradients){
          ice_gradients_ = ice_Gradients;
          assert(ice_gradients_.first.size() > 0);
          assert(ice_gradients_.second.size() == 2);
          ice_gradients_loaded_ = true;}

        template<typename DataType>
            std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const{
                assert(params.size()==27);
                //unpack things so we have legible names
                DataType convNorm=params[0];
                DataType promptNorm=params[1];
                DataType CRDeltaGamma=params[2];
                DataType piKRatio=params[3];
                DataType NeutrinoAntineutrinoRatio=params[4];
                DataType AtmosphericZenithVariationCorrectionFactorParameter=params[5];
                // barr parameters
                DataType BarrWP=params[6];
                DataType BarrWM=params[7];
                DataType BarrYP=params[8];
                DataType BarrYM=params[9];
                DataType BarrZP=params[10];
                DataType BarrZM=params[11];
                DataType BarrHP=params[12];
                DataType BarrHM=params[13];
                // detector systematics
                DataType deltaDomEff=params[14];
                DataType hqdeltaDomEff=params[15];
                DataType holeiceForward=params[16];
                DataType icegrad0=params[17];
                DataType icegrad1=params[18];
                DataType icegrad2=params[19];
                // astrophysical neutrino parameters
                DataType astroNorm=params[20];
                DataType astroDeltaGamma=params[21];
                // cross section uncertainty nuisance parameters
                DataType nuxs=params[22];
                DataType nubarxs=params[23];
                // kaon losses
                DataType kaonLosses=params[24];
                // astrophysical neutrino parameters
                DataType astroNormSec=params[25];
                DataType astroDeltaGammaSec=params[26];

                using cachedWeighter=cachedValueWeighter<DataType,Event,double>;
                cachedWeighter astroMuFlux(&Event::cachedAstroMuWeight);
                cachedWeighter convPionFlux(&Event::cachedConvPionWeight);
                cachedWeighter convKaonFlux(&Event::cachedConvKaonWeight);
                cachedWeighter convFlux(&Event::cachedConvWeight);
                cachedWeighter promptFlux(&Event::cachedPromptWeight);
                cachedWeighter BarrWPComp(&Event::cachedBarrModWP);
                cachedWeighter BarrWMComp(&Event::cachedBarrModWM);
                cachedWeighter BarrYPComp(&Event::cachedBarrModYP);
                cachedWeighter BarrYMComp(&Event::cachedBarrModYM);
                cachedWeighter BarrZPComp(&Event::cachedBarrModZP);
                cachedWeighter BarrZMComp(&Event::cachedBarrModZM);
                cachedValueWeighter<DataType,Event, float> adHocFlux(&Event::oneWeight);

                using neuaneu_t = antiparticleWeighter<Event,DataType>;
                neuaneu_t neuaneu_w(NeutrinoAntineutrinoRatio);

                using domEffW_t = DOMEffWeighter<Event,DataType>;
                domEffW_t convDOMEff(*domefficiencySplines_,deltaDomEff,FluxComponent::atmConv);
                domEffW_t promptDOMEff(*domefficiencySplines_,deltaDomEff,FluxComponent::atmPrompt);
                domEffW_t astroNuMuDOMEff(*domefficiencySplines_,deltaDomEff,FluxComponent::diffuseAstro_mu);

                using hqdomEffW_t = DOMEffWeighter<Event,DataType>;
                hqdomEffW_t hqconvDOMEff(*hqdomefficiencySplines_,hqdeltaDomEff,FluxComponent::atmConv);
                hqdomEffW_t hqpromptDOMEff(*hqdomefficiencySplines_,hqdeltaDomEff,FluxComponent::atmPrompt);

                using holeIceW_t = holeiceWeighter<Event,DataType>;
                holeIceW_t convHoleIceWeighter(*holeiceSplines_,holeiceForward,FluxComponent::atmConv);
                holeIceW_t promptHoleIceWeighter(*holeiceSplines_,holeiceForward,FluxComponent::atmPrompt);
                holeIceW_t astroNuMuHoleIceWeighter(*holeiceSplines_,holeiceForward,FluxComponent::diffuseAstro_mu);

                atmosphericDensityUncertaintyWeighter<Event,DataType> aduw(atmosphericDensityUncertaintySpline_,AtmosphericZenithVariationCorrectionFactorParameter);
                kaonLossesUncertaintyWeighter<Event,DataType> kluw(kaonLossesUncertaintySpline_,kaonLosses);

                icegradientWeighter<Event,DataType> ice_grad_weighter_0,ice_grad_weighter_1;
                if(ice_gradients_loaded_){
                  ice_grad_weighter_0 = icegradientWeighter<Event,DataType>(ice_gradients_.first,ice_gradients_.second[0],icegrad0);
                  ice_grad_weighter_1 = icegradientWeighter<Event,DataType>(ice_gradients_.first,ice_gradients_.second[1],icegrad1);
                }

                /*
                attenuationWeighter<Event,DataType> conv_nu_att_weighter(FluxComponent::atmConv, nuxs, nufate_xs_path_, true);
                attenuationWeighter<Event,DataType> conv_nubar_att_weighter(FluxComponent::atmConv, nubarxs, nufate_xs_path_, true);
                attenuationWeighter<Event,DataType> prompt_nu_att_weighter(FluxComponent::atmPrompt, nuxs, nufate_xs_path_, true);
                attenuationWeighter<Event,DataType> prompt_nubar_att_weighter(FluxComponent::atmPrompt, nubarxs, nufate_xs_path_, true);
                attenuationWeighter<Event,DataType> astro_nu_att_weighter(FluxComponent::diffuseAstro_mu, nuxs, nufate_xs_path_, true);
                attenuationWeighter<Event,DataType> astro_nubar_att_weighter(FluxComponent::diffuseAstro_mu, nubarxs, nufate_xs_path_, true);
                */

                using attW_t = attenuationWeighter<Event,DataType>;
                attW_t  conv_nu_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::atmConv), nuxs, nubarxs);
                //attW_t  conv_nubar_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::atmConv), nubarxs);
                attW_t  prompt_nu_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::atmPrompt), nuxs, nubarxs);
                //attW_t  prompt_nubar_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::atmPrompt), nubarxs);
                attW_t  astro_nu_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::diffuseAstro_mu), nuxs, nubarxs);
                //attW_t  astro_nubar_att_weighter((*attenuationSplines_), (golemfit::FluxComponent::diffuseAstro_mu), nubarxs);

                auto conventionalComponent = convNorm
                    *aduw*kluw*(convPionFlux + piKRatio*convKaonFlux + BarrWP*BarrWPComp + BarrYP*BarrYPComp + BarrZP*BarrZPComp + BarrWM*BarrWMComp + BarrYM*BarrYMComp + BarrZM*BarrZMComp)
                    *powerlawTiltWeighter<Event,DataType>(medianConvEnergy, CRDeltaGamma)
                    *convHoleIceWeighter*convDOMEff
                    *ice_grad_weighter_0
                    *ice_grad_weighter_1
                    *conv_nu_att_weighter
                ;

                auto promptComponent = promptNorm*promptFlux
                    *promptHoleIceWeighter*promptDOMEff
                    *powerlawTiltWeighter<Event,DataType>(medianPromptEnergy, CRDeltaGamma)
                    *ice_grad_weighter_0
                    *ice_grad_weighter_1
                    *prompt_nu_att_weighter
                ;

                auto astroComponent = astroNorm*astroMuFlux
                    *astroNuMuHoleIceWeighter*astroNuMuDOMEff
                    *powerlawTiltWeighter<Event,DataType>(astroPivotEnergy, astroDeltaGamma)
                    *ice_grad_weighter_0
                    *ice_grad_weighter_1
                    *astro_nu_att_weighter
                    *neuaneu_w
                ;

                auto astroLogParaboloidComponent = astroNorm*astroMuFlux
                    *logParaboloidWeighter<Event,DataType>(1e5,astroDeltaGamma,astroDeltaGammaSec)
                    *neuaneu_w
                ;

                auto adhocAstroModelComponent = astroNormSec*(adHocFlux)*
                    astroModelWeighter<Event,DataType>(steering->ad_hoc_astro_model,*astrophysicalNeutrinoModels_,astroDeltaGammaSec);

                if(fit_type == golemfit::DiffuseFitType::SinglePowerLaw) {
                    return (conventionalComponent+promptComponent+astroComponent);
                } else if(fit_type == golemfit::DiffuseFitType::LogParaboloid) {
                    return (conventionalComponent+promptComponent+astroLogParaboloidComponent);
                } else if(fit_type == golemfit::DiffuseFitType::AdHocModel){
                  if(not astrophysical_models_splines_loaded_)
                    throw std::runtime_error("Astrophysical model splines not loaded");
                  return (conventionalComponent+promptComponent+astroComponent+adhocAstroModelComponent);
                } else
                    throw std::runtime_error("Diffuse Fit type not found");
            }
};

} // close sterile namespace

#endif
