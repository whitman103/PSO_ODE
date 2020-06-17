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
#include <mpi.h>
	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/poisson_distribution.hpp>
	using boost::poisson_distribution;
#include <boost/random/variate_generator.hpp>
	using boost::variate_generator;
	
	
	
double testPass(boost::mt19937* inGenerator);
	
int main(int argc, char** argv){
	generator.seed(time(NULL));
	
	cout<<testPass(&generator)<<endl;
	cout<<(double)generator.max()<<endl;
	
	
	return 0;
}

double testPass(boost::mt19937* inGenerator){
	return (double)(*inGenerator)();
}