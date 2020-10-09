#include "utils.h"

std::istream& operator>>(std::istream& is, paramFixSpec::singleParam& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	
	is >> p.first >> c >> p.second;
	if(c!=','){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	
	is >> c;
	if(c!='}'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	return(is);
}

std::istream& operator>>(std::istream& is, paramFixSpec& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	//bool done=false;
	while(true){
		c=is.peek();
		if(is.eof())
			break;
		switch(c){
			case '}':
				is >> c;
				//done=true;
				break;
			case ',':
				is >> c;
				//fall through
			case '{':
			{
				paramFixSpec::singleParam pa;
				is >> pa;
				if(is.fail())
					return(is);
				p.params.push_back(pa);
			}
				break;
			default:
				is.setstate(is.rdstate()|std::ios::failbit);
				return(is);
		}
	}
	return(is);
}

std::vector<double> edgesToCenters(const std::vector<double> && edges) {
    std::vector<double> res(edges.size()-1);
    for(unsigned int i=0; i<edges.size()-1; ++i) {
        res[i] = (edges[i] + edges[i+1])/2.0;
    }
    return res;
}

std::vector<double> edgesToCenters(const std::vector<double> & edges) {
    std::vector<double> res(edges.size()-1);
    for(unsigned int i=0; i<edges.size()-1; ++i) {
        res[i] = (edges[i] + edges[i+1])/2.0;
    }
    return res;
}

