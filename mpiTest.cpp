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
	
	
	
	
	
int main(int argc, char** argv){
	generator.seed(time(NULL));
	
	
	
	MPI_Init(NULL,NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	cout<<generator()<<endl;sssss
	
	MPI_Finalize();
	
	return 0;
}