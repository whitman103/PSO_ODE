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
#include <filesystem>
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

	
/* Quick functions for portability. */
	
inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline int randSite(int size){
	return generator()%size;
}


#include "SIR.h"


int main(){

    string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);
    outputFolder+="\\";

    //T, I, V
    const int numOfSpecies(3);
    vector<double> speciesVector(numOfSpecies,0);
    vector<double> resetSpecies=speciesVector;
    //beta, delta, c, p
    const int numOfParameters(4);
    vector<double> initBounds={2,2,2,2};
    const int numOfParticles(20);


    vector<double (*)(Particle*,vector<double>&)> interactionFuncts(numOfSpecies);
    interactionFuncts[0]=firstInteraction;
    interactionFuncts[1]=secondInteraction;
    interactionFuncts[2]=thirdInteraction;


    double inDelta(10);
    vector<FuzzyTree> FuzzyStructure(numOfParticles,FuzzyTree(inDelta));

    vector<Particle> particleList(numOfParticles,Particle(numOfParameters,initBounds,interactionFuncts));

    vector<double> stoppingTimes={0,2,4,8,16,32};
    const int numberOfReports(5);
    double timeIncrement(0.0002);

    vector<double> fitnessHolder(numOfParticles,0);
    vector<double> globalBestSolutions(numOfParameters,0);

    for(int particle=0;particle<(int)particleList.size();particle++){
        Particle* interParticle=&particleList[particle];
        vector<vector<double> > testArray(numberOfReports,vector<double> (numOfSpecies,0));

        double currentTime(0);
        speciesVector=resetSpecies;


        for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
            rungeKuttaUpdate(interParticle,speciesVector,stoppingTimes[reportIndex],stoppingTimes[reportIndex+1],timeIncrement);
            for(int i=0;i<numOfSpecies;i++){
                testArray[reportIndex][i]=speciesVector[i];
            }
        }
        (*interParticle).currentFitness=fitnessFunction(trueArray,testArray);
        if((*interParticle).currentFitness<(*interParticle).bestFitness){
            (*interParticle).bestFitness=(*interParticle).currentFitness;
            (*interParticle).bestSolution=(*interParticle).currentSolution;
        }
        fitnessHolder[particle]=(*interParticle).currentFitness;
    }

    int bestParticle=min_element(fitnessHolder.begin(),fitnessHolder.end())-fitnessHolder.begin();
    globalBestSolutions=particleList[bestParticle].currentSolution;
    
    for(int particle=0;particle<(int)particleList.size();particle++){
        particleList[particle].performUpdate(&generator,globalBestSolutions,&FuzzyStructure[particle]);
    }


    return 0;
}