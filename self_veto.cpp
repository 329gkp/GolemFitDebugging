#include "self_veto.h"

namespace self_veto{
	
	double overburden(double cos_theta, double depth, double elevation){
		double r = 6371315 + elevation; //curvature radius of the surface (meters)
		//this is secretly a translation in polar coordinates
		double d=(cos_theta*(r-depth));
		return (sqrt(2*r*depth + d*d - depth*depth) - (r-depth)*cos_theta);
	}
	
	double effective_costheta(double costheta){
		const static double p[5]={0.102573*0.102573, -0.068287, 0.958633, 0.0407253, 0.817285};
		return sqrt((costheta*costheta + p[0] + p[1]*pow(costheta,p[2]) + p[3]*pow(costheta,p[4]))
		            /(1 + p[0] + p[1] + p[3]));
	}
	
	std::ostream& operator<<(std::ostream& os, FluxType& ft){
		switch(ft){
			case mu:
				return(os << "mu");
			case nue:
				return(os << "nue");
			case numu:
				return(os << "numu");
			case charm:
				return(os << "charm");
            default:
                return os;
		}
	}
	
	double gaisser_flux(double energy, ParticleType ptype){
		//normalizations for each element
		static const double norm[3][5]={
			{7860., 3550., 2200., 1430., 2120.},
			{20., 20., 13.4, 13.4, 13.4},
			{1.7, 1.7, 1.14, 1.14, 1.14}
		};
		//spectral indices
		static const double gamma[3][5]={
			{2.66, 2.58, 2.63, 2.67, 2.63},
			{2.4, 2.4, 2.4, 2.4, 2.4},
			{2.4, 2.4, 2.4, 2.4, 2.4},
		};
		//cutoff rigidity
		static const double rigidity[3]={4e6, 30e6, 2e9};
		
		unsigned int z=((unsigned)ptype < 100 ? 1 : (unsigned)ptype % 100);
		unsigned int idx=0;
		switch(ptype){
			case PPlus: idx=0; break;
			case He4Nucleus: idx=1; break;
			case N14Nucleus: idx=2; break;
			case Al27Nucleus: idx=3; break;
			case Fe56Nucleus: idx=4; break;
		}
		double sum=0;
		for(unsigned int i=0; i<3; i++)
			sum+=norm[i][idx]*pow(energy,(-gamma[i][idx]))*exp(-energy/(rigidity[i]*z));
		return(sum);
	}
	
	marray<double,1> logspace(double log_start, double log_stop, unsigned num){
		assert(num>0);
		auto step=(log_stop-log_start)/(num-1);
		marray<double,1> result({num});
		for(unsigned i=0; i<num; i++)
			result[i]=i;
		result*=step;
		result+=log_start;
		for(double& v : result)
			v=pow(10.,v);
		return(result);
	}
	
} //namespace self_veto
