#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include "fuzzyDef.h"
using namespace std;
#include <boost/random/mersenne_twister.hpp>
using boost::mt19937;
boost::mt19937 generator;



 
 
 
int main(){
	
	FuzzyTree testTree(3);
	
	testTree.delta=1.5;
	testTree.phi=-.75;
	testTree.setParameters();
	cout<<testTree.inertia<<endl;
	
	
	
	
	return 0;
}