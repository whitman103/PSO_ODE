#include <random>
#include <windows.h>
#include <chrono>
#include <ctime>
#include <map>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/poisson_distribution.hpp>
	using boost::poisson_distribution;
#include <boost/random/variate_generator.hpp>
	using boost::variate_generator;
	
	
/* Quick functions for portability. */
	
inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline double twoSidedRandPull(){
	return (2.*randPull()-1);
}

inline int randSite(int size){
	return generator()%size;
}



int main(){
	
	generator.seed(time(NULL));
	const int numOfSpecies(5);
	
	ofstream outStruct("randomHillStructure.txt");
	
	outStruct<<numOfSpecies<<endl;
	
	for(int i=0;i<numOfSpecies;i++){
		vector<int> interactions;
		int numOfInteractions(1+randSite(numOfSpecies-1));
		for(int j=0;j<numOfInteractions;j++){
			int interactionLabel(randSite(numOfSpecies));
			vector<int>::iterator it; 
			it=find(interactions.begin(),interactions.end(),interactionLabel);
			if(it==interactions.end()){
				interactions.push_back(interactionLabel);
			}
			else{
				j--;
			}
		}
		
		outStruct<<interactions.size()<<" ";
		for(int i=0;i<(int)interactions.size();i++){
			outStruct<<interactions[i]<<" "<<2<<" ";
		}
		outStruct<<endl;
	}
	outStruct.close();
	
	
	
	return 0;
}