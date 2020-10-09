#include <cmath>
#include <iomanip>
#include <map>
#include <vector>

#include <photospline/splinetable.h>
#include <nuSQuIDS/marray.h>
#include <boost/core/ignore_unused.hpp>

namespace self_veto{
  using namespace nusquids;
	
	///\brief Overburden for a detector buried underneath a flat surface.
	///
	///\param cos_theta cosine of zenith angle (in detector-centered coordinates)
	///\param depth depth of detector (in meters below the surface)
	///\param elevation elevation of the surface above sea level (meters)
	double overburden(double cos_theta, double depth=1950, double elevation=2400);
	
	///Minimum muon energy required to survive the given thickness of ice with
	///energy atleast emin 50% of the time.
	///
	///\param distance ice thickness (meters)
	///\param emin minimum allowed energy for surviving muons (GeV)
	///\return minimum muon energy at the surface (GeV)
	template<typename Number=double>
	Number minimum_muon_energy(double distance, Number emin=1e3){
		//coefficients must be ordered from highest degree to lowest
		auto polynomial=[](Number x, const std::vector<Number>& coefficients){
			Number sum=0;
			for(Number c : coefficients)
				sum=sum*x+c;
			return(sum);
		};
		Number le=log10(emin);
		Number a=polynomial(le,{0.187, -0.476, 2.793});
		Number b=polynomial(le,{0.023, -0.201, 2.069});
		Number c=polynomial(le,{3.882, -2.689});
		return(pow(10.,polynomial(distance,{c/1e10,b/1e4,a})));
	}
	
	///Effective local atmospheric density correction from D. Chirkin.
	///"Fluxes of atmospheric leptons at 600-GeV - 60-TeV".
	////2004. http://arxiv.org/abs/hep-ph/0407078
	double effective_costheta(double costheta);
	
	enum FluxType{
		mu=0,
		numu=1,
		nue=2,
		charm=3
	};
	
	std::ostream& operator<<(std::ostream& os, FluxType& ft);
	
	///Evaluate the lepton yield of the given lepton family per shower for showers
	///of a given type.
	///
	///\param emin lepton energy at which to evaluate
	///\param primary_energy total primary energy
	///\param primary_mass primary atomic number
	///\param cos_theta cosine of primary zenith angle
	///\param kind family of parameterization
	///\param differential if true, evaluate dN/dE at emin. Otherwise, evaluate N(>emin).
	template<typename Number=double>
	Number elbert_yield(Number emin, double primary_energy, unsigned primary_mass,
						double cos_theta, FluxType kind=mu, bool differential=false){
		struct params_type{ double a, p1, p2, p3; };
		const static std::map<FluxType,params_type> elbert_params{
			//{elbert,{14.5,0.757+1,5.25,1}},
			{mu,{49.41898933142626,0.62585930096346309+1,4.9382653076505525,0.58038589096897897}},
			{numu,{79.918830537201231,0.46284463423687988+1,4.3799061061862314,0.31657956163506323}},
			{nue,{0.548,0.669+1,8.05,0.722}},
			{charm,{780.35285355003532/1e6,-0.39555243513109928+1,7.3461490462825703,0.76688386541155051}}
		};
		const params_type& params=elbert_params.find(kind)->second;
		double En = primary_energy/primary_mass;
		Number x = emin/En;
		double decay_prob;
		if(kind == charm)
			decay_prob=1;
		else
			decay_prob = 1./(En*effective_costheta(cos_theta));
		Number icdf=0;
		if(x<1){
			icdf=params.a*primary_mass*decay_prob * pow(x,(-params.p1)) * pow(1-pow(x,params.p3),params.p2);
			if(differential)
				icdf*=(params.p1/x + params.p2*params.p3*pow(x,(params.p3-1))/(1-pow(x,params.p3)))/En;
		}
		return(icdf);
	}
	
	enum ParticleType{
		PPlus       =   14,
		He4Nucleus  =  402,
		N14Nucleus  = 1407,
		Al27Nucleus = 2713,
		Fe56Nucleus = 5626
	};
	
	///\brief Evaluate the Gaisser H3a parameterization of the cosmic ray flux.
	///
	///T. K. Gaisser. "Spectrum of cosmic-ray nucleons, kaon production, and the
	///atmospheric muon charge ratio". Astroparticle Physics, 35(12):801--806,
	///2012. ISSN 0927-6505. doi: 10.1016/j.astropartphys.2012.02.010.
	///
	///\param energy total energy of primary nucleus [GeV]
	///\param ptype: particle type code
	///\returns particle flux [particles/(GeV m^2 sr s)]
	double gaisser_flux(double energy, ParticleType ptype);
	
	//TODO: multidimensional insanity?
	marray<double,1> logspace(double log_start, double log_stop, unsigned num);
	
	//bleh only 1D case for now
	template<typename T>
	marray<T,1> diff(const marray<T,1>& a){
		marray<T,1> d({a.extent(0)-1});
		auto it=a.begin(), end=a.end();
		T prev=*it;
		it++;
		for(auto oit=d.begin(); it!=end; it++,oit++){
			*oit=*it-prev;
			prev=*it;
		}
		return(d);
	}
	
	template<typename Number=double>
	struct response{
		//contributions to the differential neutrino flux from chunks of the
		//primary spectrum for each element:
		marray<double,2> response;
		//mean integral muon yield from same chunks
		marray<Number,2> muonyield;
		marray<double,1> energy_per_nucleon;
	};
	
	///Evaluate the response function (contributions from each chunk of the cosmic
	///ray flux) for leptons of the given energy and zenith angle
	///
	///\param enu lepton energy [GeV]
	///\param emu minimum muon energy needed for veto [GeV at surface]
	///\param cos_theta: cosine of zenith angle
	///\param kind lepton type for which to evaluate the passing rate:
	///            numu for muon neutrinos from pion/kaon decay
	///            nue for electron neutrinos from kaon decay
	///            charm for either flavor from charmed meson decay
	///\return a tuple (response, muonyield, energy_per_nucleon)
	template<typename Number=double>
	response<Number> response_function(double enu, Number emu, double cos_theta, FluxType kind=numu){
		//This constant was oriignally 100, but with 1/10th the calculation we
		//still get a very good approximation
		const unsigned int nBins=10;
		response<Number> result={
			marray<double,2>({5,nBins}),
			marray<Number,2>({5,nBins}),
			marray<double,1>()
		};
		result.energy_per_nucleon = logspace(log10(enu), log10(enu)+3, nBins+1);
		marray<ParticleType,1> ptypes({5},{PPlus,He4Nucleus,N14Nucleus,Al27Nucleus,Fe56Nucleus});
		marray<int,1> A({5},{1,4,14,27,56});
		for(unsigned i=0; i<5; i++){
			marray<ParticleType,1> ptype({1},{ptypes[i]});
			auto a=A[i];
			//primary energies that contribute to the neutrino flux at given energy
			auto penergy = a*result.energy_per_nucleon;
			//width of energy bins
			auto de = diff(penergy);
			//center of energy bins
			penergy.erase(0,penergy.extent(0)-1,1);
			auto pe = penergy+de/2;
			//hobo-integrate the flux
			auto weights = applyOp(pe,ptype,[](double e,ParticleType p){return(gaisser_flux(e,p));})*de;
			marray<double,1> ey(pe.get_extents());
			for(unsigned int j=0; j<pe.extent(0); j++){
				ey[j]=elbert_yield(enu,pe[j],a,cos_theta,kind,true);
				result.response[i][j]=ey[j]*weights[j];
				result.muonyield[i][j]=elbert_yield(emu, pe[j], a, cos_theta, mu, false);
			}
		}
		//shift energy per nucleon to be at bin centers
		{
			auto de = diff(result.energy_per_nucleon);
			result.energy_per_nucleon.erase(0,result.energy_per_nucleon.extent(0)-1,1); //delete last entry
			result.energy_per_nucleon+=de/2;
		}
		return(result);
	}
	
	///Calculate the probability that neutrinos of the given energy and type
	///will be accompanied by at least one muon from an unrelated branch of the
	///shower.
	///
	///\param enu neutrino energy
	///\param emu minimum muon energy needed for veto [GeV at surface]
	///\param cos_theta cosine of zenith angle
	///\param kind neutrino type for which to evaluate the passing rate:
	///            numu for muon neutrinos from pion/kaon decay, nue for
	///            electron neutrinos from kaon decay, or charm for either
	///            flavor from charmed meson decay
	template<typename Number=double>
	Number uncorrelated_passing_rate(double enu, Number emu, double cos_theta, FluxType kind){
		//get contributions to the differential neutrino flux from chunks of
		//the cosmic ray spectrum in each element
		auto response=response_function(enu, emu, cos_theta, kind);
		auto& contrib=response.response;
		auto& muyield=response.muonyield;
		//normalize contributions over all primaries
		auto norm_contrib = contrib / contrib.sum_same_rank(contrib.rank()-1).sum_same_rank(contrib.rank()-2);
		//weight contributions by probability of having zero muons. if that
		//probability is always 1, then this returns 1 by construction
		return((exp(-muyield)*norm_contrib).sum());
	}
	
	///Calculate the differential flux of muon-neutrinos from pion and kaon decay
	///in air showers, optionally restricting the solution to the part of the
	///decay phase space where the lab-frame energy of the muon is above the given
	///threshold [Schoenert]. See [Lipari] for a detailed derivation of the
	///cascade-equation solution.
	///
	///\param enu neutrino energy at which to evaluate [GeV]
	///\param cos_theta cosine of zenith angle in detector-centered coordinates
	///\param emu if nonzero, calculate only the fraction of the flux
	///           where the partner muon has at least the given energy.
	///           Otherwise, calculate the total flux.
	///
	///\return the neutrino flux as a fraction of the primary nucleon flux
	///        evaluated at the given energy.
	///
	///[Schoenert] S. Schoenert, T. K. Gaisser, E. Resconi, and O. Schulz.
	///"Vetoing atmospheric neutrinos in a high energy neutrino telescope".
	///Phys. Rev. D, 79:043009, Feb 2009. doi: 10.1103/PhysRevD.79.043009.
	///[Lipari] P. Lipari. "Lepton spectra in the earth's atmosphere".
	///Astroparticle Physics, 1 (2):195--227, 1993. ISSN 0927-6505.
	///doi: 10.1016/0927-6505(93)90022-6.
	template<typename Number=double>
	Number analytic_numu_flux(double enu, double cos_theta, Number emu=0){
		//Spectral index of the integral nucleon flux at the top of the atmosphere
		const double GAMMA=1.7;
		//Spectrum-weighted moments for nucleon and meson production
		const double Z_NN=0.298, Z_NPI=0.079, Z_NK=0.0118, Z_PIPI=0.271, Z_KK=0.223;
        boost::ignore_unused(Z_NN);
        boost::ignore_unused(Z_NPI);
        boost::ignore_unused(Z_NK);
        boost::ignore_unused(Z_PIPI);
        boost::ignore_unused(Z_KK);
		//Critical energies for pions and kaons above which re-interaction is more
		//likely than decay in flight (GeV)
		const double EPS_PI=115., EPS_K=850.;
		const double R_PI=.5731, R_K=0.0458, ALAM_N=120., ALAM_PI=160., ALAM_K=180.;
		double F = (GAMMA + 2.)/(GAMMA + 1.);
		Number B_PI = F * (ALAM_PI - ALAM_N) / (ALAM_PI * log(ALAM_PI/ALAM_N)*(1.-R_PI));
		Number B_K = F * (ALAM_K - ALAM_N) / ( ALAM_K * log(ALAM_K /ALAM_N)*(1.-R_K));
		Number A_PI, A_K;
		if(emu!=0){
			Number z = 1 + emu/enu;
			double zpimin = 1./(1.-R_PI);
			double zkmin = 1./(1.-R_K);
			Number zzpi = (z >= zpimin ? z : zpimin);
			Number zzk = (z >= zkmin ? z : zkmin);
			A_PI = Z_NPI/((1.-R_PI)*(GAMMA+1)*pow(zzpi,(GAMMA+1.)));
			A_K =  Z_NK/((1.-R_K)*(GAMMA+1)*pow(zzk,(GAMMA+1.)));
			B_PI = zzpi*B_PI*(1.-R_PI);
			B_K = zzk*B_K*(1.-R_K);
		}
		else{
			A_PI = (Z_NPI*(pow((1.-R_PI),(GAMMA+1.)))/(GAMMA+1.))/(1.-R_PI);
			A_K =  (Z_NK *(pow((1.-R_K ),(GAMMA+1.)))/(GAMMA+1.))/(1.-R_K);
		}
		double cs=effective_costheta(cos_theta);
		return (A_PI / (1. + B_PI*cs*enu/EPS_PI) + 0.635 * A_K / (1. + B_K*cs*enu/EPS_K));
	}
	
	///Calculate the probability that muon neutrinos of the given energy will be
	///accompanied by a muon from the same decay vertex.
	///
	///\param enu neutrino energy
	///\param emu minimum muon energy needed for veto [GeV at surface]
	///\param cos_theta cosine of zenith angle
	template<typename Number=double>
	Number correlated_passing_rate(double enu, Number emu, double cos_theta){
		double flux = analytic_numu_flux(enu, cos_theta);
		Number sflux = analytic_numu_flux(enu, cos_theta, emu);
		return (flux-sflux)/flux;
	}
	
	///A combination of the Schoenert et al calculation and an approximate
	///treatment of uncorrelated muons from the rest of the shower.
	template<typename Number=double>
	class AnalyticPassingFraction{
	public:
		enum fluxKind{
			conventional,
			charm
		};
		
		static std::string cachedDataName(fluxKind kind, FluxType flux, double threshold){
			if(kind==charm)
				flux=FluxType::charm;
			std::ostringstream ss;
			ss << "jvs_data/uncorrelated_veto_prob." << flux << "."
			   << std::scientific << std::setprecision(1) << threshold << ".fits";
			return(ss.str());
		}
		
		///\param kind either conventional for neutrinos from pion/kaon decay
		///            or charm for neutrinos from charmed meson decay
		///\param veto_threshold energy at depth where a single muon is guaranteed
		///                      to be rejetected by veto cuts [GeV]
		///\param floor minimum passing fraction (helpful for avoiding numerical
		///             problems in likelihood functions)
		AnalyticPassingFraction(fluxKind kind=conventional, Number veto_threshold=1e3, double floor=1e-4):
		kind(kind),veto_threshold(veto_threshold),floor(floor),ct_min(.05)
		//nueTable(cachedDataName(kind,nue,veto_threshold)),
		//numuTable(cachedDataName(kind,numu,veto_threshold))
		{}
		
		///Estimate the fraction of atmospheric neutrinos that will arrive without
		///accompanying muons from the same air shower.
		///
		///\param particleType: neutrino type for which to evaluate the veto
		///\param enu neutrino energy [GeV]
		///\param ct cosine of the zenith angle
		///\param depth vertical depth [m]
		Number operator()(FluxType particleType, double enu, double ct, double depth) const{
			//std::cout << "log(enu) = " << log10(enu) << std::endl;
			//std::cout << "ct = " << ct << std::endl;
			//std::cout << "depth = " << depth << std::endl;
			if(ct<=ct_min)
				return(1);
			//std::cout << "mme(" << overburden(ct, depth) << ',' << veto_threshold <<')' << std::endl;
			Number emu = minimum_muon_energy(overburden(ct, depth), veto_threshold);
			//std::cout << "emu = " << emu << std::endl;
			Number pr=0;
			
			//spline implementation
			/*if(particleType==numu){
				pr=numuTable(log10(enu), ct, depth);
				//For NuMu specifically there is a guaranteed accompanying muon.
				//Estimate the passing fraction from the fraction of the decay phase
				//space where the muon has too little energy to make it to depth.
				//NB: strictly speaking this calculation applies only to 2-body
				//decays of pions and kaons, but is at least a conservative estimate
				//for 3-body decays of D mesons.
				double direct = correlated_passing_rate(enu, emu, ct);
				pr*=direct;
			}
			else if(particleType==nue){
				pr=nueTable(log10(enu), ct, depth);
			}*/
			
			//direct implementation
			if(kind==conventional){
				if(particleType==numu)
					pr=uncorrelated_passing_rate(enu, emu, ct, numu);
				else if(particleType==nue)
					pr=uncorrelated_passing_rate(enu, emu, ct, nue);
				else
					throw std::logic_error("Unexpected flux/particle combination");
			}
			else //kind==charm
				pr=uncorrelated_passing_rate(enu, emu, ct, self_veto::charm);
			
			//For NuMu specifically there is a guaranteed accompanying muon.
			//Estimate the passing fraction from the fraction of the decay phase
			//space where the muon has too little energy to make it to depth.
			//NB: strictly speaking this calculation applies only to 2-body
			//decays of pions and kaons, but is at least a conservative estimate
			//for 3-body decays of D mesons.
			if(particleType==numu){
				Number direct=correlated_passing_rate(enu, emu, ct);
				pr*=direct;
			}
			
			//std::cout << "pr = " << pr << std::endl;
			if(pr>1)
				return(1);
			if(pr<floor)
				return(floor);
			return(pr);
		}
		
	private:
		fluxKind kind;
		Number veto_threshold;
		double floor;
		double ct_min;
		//Splinetable nueTable;
		//Splinetable numuTable;
    public:
        double get_ct_min() {
            return ct_min;
        }
        double get_floor() {
            return floor;
        }
	};
	
} //namespace self_veto


struct ThresholdCorrection{
	ThresholdCorrection():
	nue_2010("jvs_data/threshold_correction/efficiency.nue.2010.fits"),
	numu_2010("jvs_data/threshold_correction/efficiency.numu.2010.fits")
	{}
	
	///Calculate the relative efficiency at a given deposited energy as a function of
	///energy scale shift from the ratio of the effective volume at *energy x shift* to
	//the effective volume at *energy*.
	///
	///\param year data-taking year
	///\param ptype particle type as a PDG code. Only NuMu (abs(ptype)==14) are
	///             treated differently; NuTau behave similarly to NuE
	///\param energy deposited energy at bin center. The bins should be in _target_
	///              DOM efficiency, i.e. energies shifted from the native value.
	///\param domeff target DOM efficiency
	///\param base_domeff native DOM efficiency of the simulation
	double operator()(int year, int ptype, double energy, double domeff, double base_domeff){
		auto useSpline=[&](const photospline::splinetable<>& spline){
			double lo=spline.lower_extent(0);
			double hi=spline.upper_extent(0);
			double loge=log10(energy);
			double dloge=log10(domeff/base_domeff);
			double x1;
			if(loge<lo) x1=lo;
			else if(loge>hi-dloge) x1=hi-dloge;
			else x1=loge;
			double x2=x1+dloge;
			double y1=spline(&x1);
			double y2=spline(&x2);
			//std::cout << " f(" << x1 << ")=" << y1 << " f(" << x2 << ")=" << y2 << std::endl;
			return(y2/y1);
		};
		if(year==2010){
			if(std::abs(ptype)==14)
				return(useSpline(numu_2010));
			else
				return(useSpline(nue_2010));
		}
		else
			throw std::runtime_error("Tables not loaded for year "+std::to_string(year));
	}
private:
	photospline::splinetable<> nue_2010, numu_2010;
};
