#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <boost/crc.hpp>
#include <boost/shared_ptr.hpp>
#include <hdf5.h>
#include <nuSQuIDS/tools.h>
#include <PhysTools/autodiff.h>

//Stuff to enable converting between types of shared_ptrs
namespace {
	template<class SharedPointer> struct Holder {
		SharedPointer p;

		Holder(const SharedPointer &p) : p(p) {}
		Holder(const Holder &other) : p(other.p) {}
		//Holder(Holder &&other) : p(std::move<SharedPointer>(other.p)) {}

		void operator () (...) const {}
    };

	template<class T> std::shared_ptr<T> to_std_ptr(const boost::shared_ptr<T> &p) {
		typedef Holder<std::shared_ptr<T>> H;
		if(H* h = boost::get_deleter<H, T>(p))
			return h->p;
		return std::shared_ptr<T>(p.get(), Holder<boost::shared_ptr<T>>(p));
	}

	template<class T> boost::shared_ptr<T> to_boost_ptr(const std::shared_ptr<T> &p){
		typedef Holder<boost::shared_ptr<T>> H;
		if(H* h = std::get_deleter<H, T>(p))
			return h->p;
		return boost::shared_ptr<T>(p.get(), Holder<std::shared_ptr<T>>(p));
	}
}

namespace golemfit {
namespace detail {
template<unsigned int angle, typename T>
struct convertFlavorRatioToAngle {
};

template<typename T>
struct convertFlavorRatioToAngle<0, T> {
    T operator()(T& fr0, T& fr1, T& fr2) {
        T a = (fr2 - 1.0)*(fr2 - 1.0);
        if(a < 0)
            a = 0;
        else if(a > 1)
            a = 1;
        return a;
    }
};

template<typename T>
struct convertFlavorRatioToAngle<1, T> {
    T operator()(T& fr0, T& fr1, T& fr2) {
        if(fr2 - 1.0 == T(0.0))
            return T(0.0);
        T b = (fr1*2.0 + fr2 - 1.0)*(fr2 - 1.0);
        if(b < -1)
            b = -1;
        else if(b > 1)
            b = 1;
        return b;
    }
};

template<unsigned int flavor_ratio, typename T>
struct convertAngleToFlavorRatio {
};

template<typename T>
struct convertAngleToFlavorRatio<0, T> {
    T operator()(T& a, T& b) {
        T res = (b + 1.0)*(sqrt(a))/2.0;
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        return res;
    }
};

template<unsigned int Dim, typename T>
struct convertAngleToFlavorRatio<0, phys_tools::autodiff::FD<Dim, T>> {
    using result_type=phys_tools::autodiff::FD<Dim, T>;
    result_type operator()(result_type& a, result_type& b) {
        unsigned int didx_a = Dim;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(a.derivative(i)!=0){
                didx_a=i;
                break;
            }
        }
        assert(didx_a < Dim);

        unsigned int didx_b = Dim;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(b.derivative(i)!=0){
                didx_b=i;
                break;
            }
        }
        assert(didx_b < Dim);

        result_type res((b.value() + 1.0)*(sqrt(a.value()))/2.0);
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        if(a.value()==0)
            res.setDerivative(didx_a, std::numeric_limits<T>::max());
        else
            res.setDerivative(didx_a, (b.value()+1.0)/(4.0*sqrt(a.value())));
        res.setDerivative(didx_b, sqrt(a.value())/2.0);
        return res;
    }
};

template<typename T>
struct convertAngleToFlavorRatio<1, T> {
    T operator()(T& a, T& b) {
        T res = (-b + 1.0)*(sqrt(a))/2.0;
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        return res;
    }
};

template<unsigned int Dim, typename T>
struct convertAngleToFlavorRatio<1, phys_tools::autodiff::FD<Dim, T>> {
    using result_type=phys_tools::autodiff::FD<Dim, T>;
    result_type operator()(result_type& a, result_type& b) {
        unsigned int didx_a;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(a.derivative(i)!=0){
                didx_a=i;
                break;
            }
        }
        unsigned int didx_b;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(b.derivative(i)!=0){
                didx_b=i;
                break;
            }
        }

        result_type res((-b.value() + 1.0)*(sqrt(a.value()))/2.0);
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        if(a.value()==0)
            res.setDerivative(didx_a, std::numeric_limits<T>::max());
        else
            res.setDerivative(didx_a, (-b.value()+1.0)/(4.0*sqrt(a.value())));
        res.setDerivative(didx_b, -sqrt(a.value())/2.0);
        return res;
    }
};

template<typename T>
struct convertAngleToFlavorRatio<2, T> {
    T operator()(T& a, T& b) {
        T res = -sqrt(a) + 1.0;
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        return res;
    }
};

template<unsigned int Dim, typename T>
struct convertAngleToFlavorRatio<2, phys_tools::autodiff::FD<Dim, T>> {
    using result_type=phys_tools::autodiff::FD<Dim, T>;
    result_type operator()(result_type& a, result_type& b) {
        unsigned int didx_a;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(a.derivative(i)!=0){
                didx_a=i;
                break;
            }
        }
        unsigned int didx_b;
        for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
            if(b.derivative(i)!=0){
                didx_b=i;
                break;
            }
        }

        result_type res(-sqrt(a.value())+1.0);
        if(res < 0)
            res = 0;
        else if(res > 1)
            res = 1;
        if(a.value()==0)
            res.setDerivative(didx_a, std::numeric_limits<T>::max());
        else
            res.setDerivative(didx_a, -1.0/(2.0*sqrt(a.value())));
        res.setDerivative(didx_b, 0);
        return res;
    }
};

} // namespace detail
} // namespace golemfit

///Wraps a vector of doubles for convenient printing
template<typename T>
struct prettyPrint{
private:
	std::vector<T>& v;
public:
	prettyPrint(std::vector<T>& v):v(v){}

	template<typename U>
	friend std::ostream& operator<<(std::ostream&, const prettyPrint<U>&);
	template<typename U>
	friend std::istream& operator>>(std::istream&, prettyPrint<U>&);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const prettyPrint<T>& p){
	bool first=true;
	os << '{';
	for(auto x : p.v){
		if(!first)
			os << ", ";
		else
			first=false;
		os << x;
	}
	return(os << '}');
}
template<typename T>
std::istream& operator>>(std::istream& is, prettyPrint<T>& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	bool done=false;
	std::vector<T> buffer;
	while(!done){
		c=is.peek();
		if(is.eof())
			break;
		switch(c){
			case '}':
				is >> c;
				done=true;
				break;
			case ',':
				is >> c;
				//fall through
			default:
			{
				T d;
				is >> d;
				if(is.fail())
					return(is);
				buffer.push_back(d);
			}
				break;
		}
	}
	p.v.swap(buffer);
	return(is);
}

///A container to help with parsing index/value pairs for fixing likelihood parameters from the commandline
struct paramFixSpec{
	struct singleParam : std::pair<unsigned int, double>{};
	std::vector<std::pair<unsigned int,double>> params;
};

std::istream& operator>>(std::istream& is, paramFixSpec::singleParam& p);
std::istream& operator>>(std::istream& is, paramFixSpec& p);

std::vector<double> edgesToCenters(const std::vector<double> && edges);
std::vector<double> edgesToCenters(const std::vector<double> & edges);

#endif //UTILS_H
