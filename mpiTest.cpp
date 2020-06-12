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
	double toSend(0);
	double receiveTest[5];
	int taskId(0);
	
	MPI_Init(NULL,NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	generator.seed(taskId*36);
	
	for(int i=0;i<world_size;i++){
		toSend=(double)generator()/(double)generator.max();
	}
	MPI_Gather(&toSend, 1, MPI_DOUBLE, receiveTest,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(taskId==0){
		for(int i=0;i<world_size;i++){
			if(toSend<receiveTest[i]){
				toSend=receiveTest[i];
			}
		}
		cout<<" master value is: "<<toSend<<endl;
	}
	
	MPI_Bcast(&toSend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	cout<<"slave value is: "<<toSend<<endl;
	
	
	MPI_Finalize();
	
	return 0;
}