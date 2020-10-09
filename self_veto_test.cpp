#include <iostream>
#include "self_veto.h"
#include <PhysTools/autodiff.h>

using namespace phys_tools::autodiff;

void test_overburden(){
	std::cout << "overburden(0.1): \n expect: " << 19213.306557239266 << "\n got:    ";
	std::cout << self_veto::overburden(0.1) << std::endl;
	std::cout << "overburden(1): \n expect: " << 1950.0 << "\n got:    ";
	std::cout << self_veto::overburden(1.0) << std::endl;
	std::cout << "overburden(0.5,2500): \n expect: " << 4997.0616867425852 << "\n got:    ";
	std::cout << self_veto::overburden(0.5,2500) << std::endl;
	std::cout << "overburden(0.5,1950,3000): \n expect: " << 3898.2120363349095 << "\n got:    ";
	std::cout << self_veto::overburden(0.5,1950,3000) << std::endl;
}

void test_minimum_muon_energy(){
	std::cout << "minimum_muon_energy(2000): \n expect: " << 2433.2279697613926 << "\n got:    ";
	std::cout << self_veto::minimum_muon_energy(2000) << std::endl;
	std::cout << "minimum_muon_energy(10000): \n expect: " << 64650.218944342829 << "\n got:    ";
	std::cout << self_veto::minimum_muon_energy(10000) << std::endl;
	std::cout << "minimum_muon_energy(2000,2e3): \n expect: " << 3932.6495801287356 << "\n got:    ";
	std::cout << self_veto::minimum_muon_energy(2000,2e3) << std::endl;
	std::cout << "minimum_muon_energy(10000,1e2): \n expect: " << 25046.670357235907 << "\n got:    ";
	std::cout << self_veto::minimum_muon_energy(10000,1e2) << std::endl;
	
	{
		FD<1> muonThresh(1e3,0);
		auto minEn=self_veto::minimum_muon_energy(5000,muonThresh);
		double eps=2;
		double minEn2=self_veto::minimum_muon_energy(5000,(double)muonThresh+eps);
		std::cout << "differentiability:\n";
		std::cout << " base value: " << minEn << '\n';
		std::cout << " extrapolated value: " << (double)minEn+eps*minEn.derivative(0) << '\n';
		std::cout << " actual other value: " << minEn2 << std::endl;
	}
}

void test_effective_costheta(){
	std::cout << "effective_costheta(0.1): \n expect: " << 0.13980645164973715 << "\n got:    ";
	std::cout << self_veto::effective_costheta(0.1) << std::endl;
	std::cout << "effective_costheta(0.5): \n expect: " << 0.50279659755690387 << "\n got:    ";
	std::cout << self_veto::effective_costheta(0.5) << std::endl;
	std::cout << "effective_costheta(0.8): \n expect: " << 0.80014429685682797 << "\n got:    ";
	std::cout << self_veto::effective_costheta(0.8) << std::endl;
}

void test_elbert_yield(){
	std::cout << "elbert_yield(1e2,1e4,1,0.1): \n expect: " << 44.32122051874078 << "\n got:    ";
	std::cout << self_veto::elbert_yield(1e2,1e4,1,0.1) << std::endl;
	std::cout << "elbert_yield(1e2,1e4,1,0.1,mu,true): \n expect: " << 0.8148351198632585 << "\n got:    ";
	std::cout << self_veto::elbert_yield(1e2,1e4,1,0.1,self_veto::mu,true) << std::endl;
	std::cout << "elbert_yield(1e2,1e4,1,0.2,numu,false): \n expect: " << 9.667977033007869 << "\n got:    ";
	std::cout << self_veto::elbert_yield(1e2,1e4,1,0.2,self_veto::numu,false) << std::endl;
	std::cout << "elbert_yield(1e3,1e5,1,0.3,nue,true): \n expect: " << 0.029092502943959972 << "\n got:    ";
	std::cout << self_veto::elbert_yield(1e2,1e5,1,0.3,self_veto::nue,true) << std::endl;
	std::cout << "elbert_yield(1e3,1e5,4,0.5,charm,false): \n expect: " << 0.011400708338409714 << "\n got:    ";
	std::cout << self_veto::elbert_yield(1e3,1e5,4,0.5,self_veto::charm,false) << std::endl;
	{
		FD<1> emin(1e2,0);
		auto yield=self_veto::elbert_yield(emin,1e4,1,0.1);
		double eps=0.5;
		double yield2=self_veto::elbert_yield((double)emin+eps,1e4,1,0.1);
		std::cout << "differentiability:\n";
		std::cout << " base value: " << yield << '\n';
		std::cout << " extrapolated value: " << (double)yield+eps*yield.derivative(0) << '\n';
		std::cout << " actual other value: " << yield2 << std::endl;
	}
}

void test_gaisser_flux(){
	std::cout << "gaisser_flux(1e3,PPlus): \n expect: " << 8.3652865693395337e-05 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e3,self_veto::PPlus) << std::endl;
	std::cout << "gaisser_flux(1e4,PPlus): \n expect: " << 1.8506172134891679e-07 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e4,self_veto::PPlus) << std::endl;
	std::cout << "gaisser_flux(1e5,PPlus): \n expect: " << 4.0584028043986309e-10 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e5,self_veto::PPlus) << std::endl;
	std::cout << "gaisser_flux(1e3,He4Nucleus): \n expect: " << 6.5960462441726051e-05 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e3,self_veto::He4Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e4,He4Nucleus): \n expect: " << 1.751513786239588e-07 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e4,self_veto::He4Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e5,He4Nucleus): \n expect: " << 4.6303346208331554e-10 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e5,self_veto::He4Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e3,N14Nucleus): \n expect: " << 2.9257885906884579e-05 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e3,self_veto::N14Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e4,N14Nucleus): \n expect: " << 7.0067336281843135e-08 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e4,self_veto::N14Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e5,N14Nucleus): \n expect: " << 1.6972643398187752e-10 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e5,self_veto::N14Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e3,Al27Nucleus): \n expect: " << 1.4891633330222113e-05 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e3,self_veto::Al27Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e4,Al27Nucleus): \n expect: " << 3.3523386115819009e-08 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e4,self_veto::Al27Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e5,Al27Nucleus): \n expect: " << 7.8289593885937939e-11 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e5,self_veto::Al27Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e3,Fe56Nucleus): \n expect: " << 2.822803878597996e-05 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e3,self_veto::Fe56Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e4,Fe56Nucleus): \n expect: " << 6.7669060376749999e-08 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e4,self_veto::Fe56Nucleus) << std::endl;
	std::cout << "gaisser_flux(1e5,Fe56Nucleus): \n expect: " << 1.6447854359180513e-10 << "\n got:    ";
	std::cout << self_veto::gaisser_flux(1e5,self_veto::Fe56Nucleus) << std::endl;
}

void test_response_function(){
	auto response=self_veto::response_function(1e4,1e3,1.0);
	std::cout << response.response.extent(0) << ',' << response.response.extent(1) << ' '
	          << response.muonyield.extent(0) << ',' << response.muonyield.extent(1) << ' '
	          << response.energy_per_nucleon.extent(0) << '\n';
	std::cout << ' ' << response.response[1][1] << ' ' << response.response[3][3] << std::endl;
	std::cout << ' ' << response.muonyield[1][1] << ' ' << response.muonyield[3][3] << std::endl;
	std::cout << ' ' << response.energy_per_nucleon[1] << ' ' << response.energy_per_nucleon[3] << std::endl;
}

void test_uncorrelated_passing_rate(){
	std::cout << "uncorrelated_passing_rate(1e5,1e3,0.8,numu): \n expect: " << 0.034798667389786278 << "\n got:    ";
	std::cout << self_veto::uncorrelated_passing_rate(1e5,1e3,0.8,self_veto::numu) << std::endl;
	std::cout << "uncorrelated_passing_rate(1e5,1e4,0.6,nue): \n expect: " << 0.74134289409280762 << "\n got:    ";
	std::cout << self_veto::uncorrelated_passing_rate(1e5,1e4,0.6,self_veto::nue) << std::endl;
	std::cout << "uncorrelated_passing_rate(1e4,1e3,0.4,charm): \n expect: " << 0.17958096874471324 << "\n got:    ";
	std::cout << self_veto::uncorrelated_passing_rate(1e4,1e3,0.4,self_veto::charm) << std::endl;
	{
		FD<1> emu(1e3,0);
		auto rate=self_veto::uncorrelated_passing_rate(1e4,emu,0.4,self_veto::charm);
		double eps=10;
		double rate2=self_veto::uncorrelated_passing_rate(1e4,(double)emu+eps,0.4,self_veto::charm);
		std::cout << "differentiability:\n";
		std::cout << " base value: " << rate << '\n';
		std::cout << " extrapolated value: " << (double)rate+eps*rate.derivative(0) << '\n';
		std::cout << " actual other value: " << rate2 << std::endl;
	}
}

void test_analytic_numu_flux(){
	std::cout << "analytic_numu_flux(1e3,0.2): \n expect: " << 0.0030596447767008367 << "\n got:    ";
	std::cout << self_veto::analytic_numu_flux(1e3,0.2) << std::endl;
	std::cout << "analytic_numu_flux(1e5,0.6): \n expect: " << 3.5030250166819709e-05 << "\n got:    ";
	std::cout << self_veto::analytic_numu_flux(1e5,0.6) << std::endl;
	std::cout << "analytic_numu_flux(1e3,0.2,1e2): \n expect: " << 0.0027989457693234541 << "\n got:    ";
	std::cout << self_veto::analytic_numu_flux(1e3,0.2,1e2) << std::endl;
	std::cout << "analytic_numu_flux(1e5,0.6,1e4): \n expect: " << 3.0071349219385648e-05 << "\n got:    ";
	std::cout << self_veto::analytic_numu_flux(1e5,0.6,1e4) << std::endl;
	{
		FD<1> emu(1e2,0);
		auto flux=self_veto::analytic_numu_flux(1e3,0.2,emu);
		double eps=1;
		double flux2=self_veto::analytic_numu_flux(1e3,0.2,(double)emu+eps);
		std::cout << "differentiability:\n";
		std::cout << " base value: " << flux << '\n';
		std::cout << " extrapolated value: " << (double)flux+eps*flux.derivative(0) << '\n';
		std::cout << " actual other value: " << flux2 << std::endl;
	}
}

void test_correlated_passing_rate(){
	std::cout << "correlated_passing_rate(1e4,1e2,1.0): \n expect: " << 0.0 << "\n got:    ";
	std::cout << self_veto::correlated_passing_rate(1e4,1e2,1.0) << std::endl;
	std::cout << "correlated_passing_rate(1e4,1e3,1.0): \n expect: " << 0.138622111420283 << "\n got:    ";
	std::cout << self_veto::correlated_passing_rate(1e4,1e3,1.0) << std::endl;
	std::cout << "correlated_passing_rate(1e3,1e2,1.0): \n expect: " << 0.11735472969218033 << "\n got:    ";
	std::cout << self_veto::correlated_passing_rate(1e3,1e2,1.0) << std::endl;
	std::cout << "correlated_passing_rate(1e4,1e3,0.5): \n expect: " << 0.13544567751558256 << "\n got:    ";
	std::cout << self_veto::correlated_passing_rate(1e4,1e3,0.5) << std::endl;
	{
		FD<1> emu(1e3,0);
		auto rate=self_veto::correlated_passing_rate(1e4,emu,0.5);
		double eps=10;
		double rate2=self_veto::correlated_passing_rate(1e4,(double)emu+eps,0.5);
		std::cout << "differentiability:\n";
		std::cout << " base value: " << rate << '\n';
		std::cout << " extrapolated value: " << (double)rate+eps*rate.derivative(0) << '\n';
		std::cout << " actual other value: " << rate2 << std::endl;
	}
}

void test_AnalyticPassingFraction(){
	{
		self_veto::AnalyticPassingFraction<> apf; //all defaults
		std::cout << "AnalyticPassingFraction()(numu,1e4,0.8,2000): \n expect: " << 0.3290739779590018 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(numu,1e5,0.8,2000): \n expect: " << 0.0001 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e5,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(numu,1e4,0.5,2000): \n expect: " << 0.5203868033516635 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(numu,1e4,0.5,1500): \n expect: " << 0.36477748900336854 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
		std::cout << "AnalyticPassingFraction()(nue,1e4,0.8,2000): \n expect: " << 0.6811762452125549 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(nue,1e5,0.8,2000): \n expect: " << 0.3305498957633972 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e5,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(nue,1e4,0.5,2000): \n expect: " << 0.7914738059043884 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction()(nue,1e4,0.5,1500): \n expect: " << 0.6681484580039978 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
	}
	{
		self_veto::AnalyticPassingFraction<> apf(self_veto::AnalyticPassingFraction<>::charm);
		std::cout << "AnalyticPassingFraction(charm)(numu,1e4,0.8,2000): \n expect: " << 0.34688885014442555 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(numu,1e5,0.8,2000): \n expect: " << 0.0001 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e5,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(numu,1e4,0.5,2000): \n expect: " << 0.5412277387847783 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(numu,1e4,0.5,1500): \n expect: " << 0.3861259980474478 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
		std::cout << "AnalyticPassingFraction(charm)(nue,1e4,0.8,2000): \n expect: " << 0.7531234622001648 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(nue,1e5,0.8,2000): \n expect: " << 0.3858462870121002 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e5,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(nue,1e4,0.5,2000): \n expect: " << 0.8493359088897705 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(charm)(nue,1e4,0.5,1500): \n expect: " << 0.7449641227722168 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
	}
	{
		self_veto::AnalyticPassingFraction<> apf(self_veto::AnalyticPassingFraction<>::conventional,1e2);
		std::cout << "AnalyticPassingFraction(conventional,1e2)(numu,1e4,0.8,2000): \n expect: " << 0.04957910352982559 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(numu,1e3,0.8,2000): \n expect: " << 0.5298359700530834 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e3,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(numu,1e4,0.5,2000): \n expect: " << 0.15561071989962025 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(numu,1e4,0.5,1500): \n expect: " << 0.06238007686187015 << "\n got:    ";
		std::cout << apf(self_veto::numu,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
		std::cout << "AnalyticPassingFraction(conventional,1e2)(nue,1e4,0.8,2000): \n expect: " << 0.27366116642951965 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(nue,1e3,0.8,2000): \n expect: " << 0.7201792001724243 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e3,0.8,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(nue,1e4,0.5,2000): \n expect: " << 0.42904481291770935 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,2000) << std::endl;
		std::cout << "AnalyticPassingFraction(conventional,1e2)(nue,1e4,0.5,1500): \n expect: " << 0.25699785351753235 << "\n got:    ";
		std::cout << apf(self_veto::nue,1e4,0.5,1500) << std::endl;
		std::cout << "--\n";
	}
	{
		FD<1> thresh(1e3,0);
		self_veto::AnalyticPassingFraction<FD<1>> apf(self_veto::AnalyticPassingFraction<FD<1>>::conventional,thresh);
		auto frac=apf(self_veto::numu,1e4,0.5,1500);
		double eps=10;
		self_veto::AnalyticPassingFraction<> apf2(self_veto::AnalyticPassingFraction<>::conventional,(double)thresh+eps);
		double frac2=apf2(self_veto::numu,1e4,0.5,1500);
		
		std::cout << "differentiability:\n";
		std::cout << " base value: " << frac << '\n';
		std::cout << " extrapolated value: " << (double)frac+eps*frac.derivative(0) << '\n';
		std::cout << " actual other value: " << frac2 << std::endl;
	}
	{
		self_veto::AnalyticPassingFraction<> apf; //all defaults
		volatile double dummy;
		std::chrono::high_resolution_clock::time_point t1, t2;
		
		unsigned count=10000;
		t1 = std::chrono::high_resolution_clock::now();
		for(unsigned i=0; i<count; i++)
			dummy=apf(self_veto::numu,1e4,0.8,2000);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()/count << " seconds/eval" << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		for(unsigned i=0; i<count; i++)
			dummy=apf(self_veto::nue,1e4,0.8,2000);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()/count << " seconds/eval" << std::endl;
	}
}

int main(){
	std::cout.precision(16);
	test_overburden();
	std::cout << "----\n";
	test_minimum_muon_energy();
	std::cout << "----\n";
	test_effective_costheta();
	std::cout << "----\n";
	test_elbert_yield();
	std::cout << "----\n";
	test_gaisser_flux();
	std::cout << "----\n";
	test_response_function();
	std::cout << "----\n";
	test_uncorrelated_passing_rate();
	std::cout << "----\n";
	test_analytic_numu_flux();
	std::cout << "----\n";
	test_correlated_passing_rate();
	std::cout << "----\n";
	test_AnalyticPassingFraction();
	
	return(0);
}