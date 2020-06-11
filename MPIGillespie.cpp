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
#include "GillespieFunctions.h"
#include <mpi.h>
using namespace std;
#include <boost/random/mersenne_twister.hpp>
  using boost::mt19937;
 boost::mt19937 generator;


inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline double twoSidedRandPull(){
	return (2.*randPull()-1);
}

inline int randSite(int size){
	return generator()%size;
}

inline double invHill(double constant, double power, double argument, double norm){
	return norm*pow(constant,power)/(pow(constant,power)+pow(argument,power));
}

inline double Hill(double constant, double power, double argument,double norm){
	return norm*pow(argument,power)/(pow(constant,power)+pow(argument,power));
}

typedef struct{
	bool hillBool;
	double constant;
	double power;
	double normalization;
	int speciesLabel;
} hillStruct;

typedef struct{
	vector<vector<hillStruct> > sampleSolution;
	vector<vector<tuple<double,double> > > bestPosition;
	vector<vector<tuple<double,double> > > currentPosition; //Normalization constant
	vector<vector<tuple<double,double> > > currentVelocity;
	vector<double> decayConsts;
	vector<double> decayVelocities;
	vector<double> bestDecayConsts;
	double bestWellness;
	double currentWellness;
} particle;

tuple<int,vector<double>> DistributionFromValues(vector<double> Values);
double returnMean(vector<double> inDist);
void loadHillStructDetails(vector<vector<hillStruct> >& inStructure, string inFile);
string convertHillStructToReactions(vector<vector<hillStruct> >& outStruct);
		

int main(int argc, char* argv[]){
	
	string baseHillStructure("randomHillStructure");
	
	vector<vector<hillStruct> > hillDetails;
	loadHillStructDetails(hillDetails, baseHillStructure);
	
	string structureLabel=convertHillStructToReactions(hillDetails);
	
	Gillespie ReactionObject1(structureLabel+"Coeffs");
	generator.seed(time(NULL));
	vector<int> specNum;
	vector<int> resetSpecNum;
	ReactionObject1.initializeData(structureLabel+"Consts",ReactionObject1.reactConsts,specNum);
	vector<double> resetConsts=ReactionObject1.reactConsts;
	resetSpecNum=specNum;
	
	
	int taskID(-1);
	int nTasks(-1);
	int numOfRuns(100);
	double gamma(1./2000.);
	
	int numParticles(2);
	vector<particle> particleSwarm;
	
	for(int i=0;i<numParticles;i++){
		particle interParticle;
		interParticle.sampleSolution=hillDetails;
	}
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	
	if(nTasks!=numParticles){
		cout<<"Reset Particle Size"<<endl;
		return 0;
	}
	generator.seed(time(NULL)+taskID);
	
	
	
	vector<double> outDist(numOfRuns,0);

	for(int i=0;i<numOfRuns;i++){
		specNum=resetSpecNum;
		double runTime(0);
		
		do{
			int reactionIndex(0);
			
			for(int i=0;i<(int)hillDetails.size();i++){
				for(int j=0;j<(int)hillDetails[i].size();j++){
					ReactionObject1.reactConsts[reactionIndex]*=Hill(hillDetails[i][j].constant,hillDetails[i][j].power, specNum[hillDetails[i][j].speciesLabel],hillDetails[i][j].normalization);
					reactionIndex++;
				}
			}
			for(int i=0;i<(int)specNum.size();i++){
				specNum[i]=currentParticle.decayConsts[i];
			}
			
			
			tuple<int,double> hold=ReactionObject1.PerformTimeStep2(specNum);
			runTime+=get<1>(hold);
			ReactionObject1.specChange(specNum,get<0>(hold),ReactionObject1.changeCoeffs);
			if(runTime>120){
				outDist[i]=specNum[10];
			}
		}while(runTime<150);
			
			//Distribution work goes here
			
	}
	
	double c1(2), c2(2);
	particle* currentParticle=&particleSwarm[sol];
	for(int i=0;i<(int)(*currentParticle).currentVelocity.size();i++){
		for(int j=0;j<(int)(*currentParticle).currentVelocity[i].size();j++){
			double rand1(randPull());
			double rand2(randPull());
			auto[currentNorm,currentConst]=(*currentParticle).currentPosition[i][j];
			auto[curBestNorm,curBestConst]=(*currentParticle).bestPosition[i][j];
			auto[groupBestNorm,groupBestConst]=particleSwarm[indexOfBestFit].currentPosition[i][j];
			double proposedNormUpdate(c1*rand1*(curBestNorm-currentNorm)+c2*rand2*(groupBestNorm-currentNorm));
			if(proposedNormUpdate+get<0>((*currentParticle).currentVelocity[i][j])>2){
				get<0>((*currentParticle).currentVelocity[i][j])=2;
			}
			else{
				get<0>((*currentParticle).currentVelocity[i][j])+=proposedNormUpdate;
			}
			if(proposedNormUpdate+get<0>((*currentParticle).currentVelocity[i][j])<-2){
				get<0>((*currentParticle).currentVelocity[i][j])=-2;
			}
			double proposedConstUpdate=c1*rand1*(curBestConst-currentConst)+c2*rand2*(groupBestConst-currentConst);
			if(proposedConstUpdate+get<1>((*currentParticle).currentVelocity[i][j])>2){
				get<1>((*currentParticle).currentVelocity[i][j])=2;
			}
			else{
				get<1>((*currentParticle).currentVelocity[i][j])+=proposedConstUpdate;
			}
			if(proposedConstUpdate+get<1>((*currentParticle).currentVelocity[i][j])<-2){
				get<1>((*currentParticle).currentVelocity[i][j])=-2;
			}
			get<0>((*currentParticle).currentPosition[i][j])+=get<0>((*currentParticle).currentVelocity[i][j]);
			get<1>((*currentParticle).currentPosition[i][j])+=get<1>((*currentParticle).currentVelocity[i][j]);
		}
	}
	
	for(int i=0;i<(int)specNum.size();i++){
		double rand1(randPull());
		double rand2(randPull());
		
		
		reactConsts[i]=currentParticle.decayConsts[i];
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Finalize();
	

	
	
	return 0;
}

string convertHillStructToReactions(vector<vector<hillStruct> >& outStruct){
	string outLeader("Structure");
	int numReactions(0);
	for(int i=0;i<(int)outStruct.size();i++){
		for(int j=0;j<(int)outStruct[i].size();j++){
			numReactions++;
		}
	}
	ofstream StructCoeffs(outLeader+"Coeffs.txt");
	StructCoeffs<<outStruct.size()<<endl;
	StructCoeffs<<(numReactions+outStruct.size())<<endl;
	//Puts out a list prop coeffs which will premultiply the propensities with the number of the concerned species, then prints the decay constants
	for(int i=0;i<(int)outStruct.size();i++){
		for(int j=0;j<(int)outStruct[i].size();j++){
			StructCoeffs<<1<<" "<<i<<" "<<1<<endl;
		}
	}
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<1<<" "<<i<<" "<<1<<endl;
	}
	
	for(int i=0;i<(int)outStruct.size();i++){
		for(int j=0;j<(int)outStruct[i].size();j++){
			StructCoeffs<<1<<" "<<i<<" "<<-1<<endl;
		}
	}
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<1<<" "<<i<<" "<<1<<endl;
	}
	
	StructCoeffs.close();
	StructCoeffs.open(outLeader+"Consts.txt");
	
	StructCoeffs<<outStruct.size()<<endl;
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<randSite(50)<<endl;
	}
	StructCoeffs<<endl;
	StructCoeffs<<(numReactions+outStruct.size())<<endl;
	StructCoeffs<<endl;
	for(int i=0;i<(int)outStruct.size();i++){
		for(int j=0;j<(int)outStruct[i].size();j++){
			StructCoeffs<<1<<endl;
		}
	}
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<1<<endl;
	}
	StructCoeffs.close();
	
	return outLeader;
	
}



double returnMean(vector<double> inList){
	double outHold(0);
	for(int i=0;i<(int)inList.size();i++){
		outHold+=inList[i]/(double)inList.size();
	}
	return outHold;
}

tuple<int,vector<double>> DistributionFromValues(vector<double> Values){
	int Counts=Values.size();
	vector<double> outputDistribution(*max_element(Values.begin(),Values.end())+1,0);
	for(int i=0;i<(int)Values.size();i++){
		outputDistribution[Values[i]]+=1/(double)Counts;
	}
	return make_tuple(Counts,outputDistribution);
}

void loadHillStructDetails(vector<vector<hillStruct> >& inStructure, string inFile){
	ifstream inData(inFile+".txt");
	if(!inData.good()){
		cout<<"Failed to open "+inFile+".txt"<<endl;
	}
	int numberOfSpecies(0);
	inData>>numberOfSpecies;
	for(int i=0;i<numberOfSpecies;i++){
		int numberOfInteractions(0);
		inData>>numberOfInteractions;
		vector<hillStruct> structureVector(numberOfInteractions);
		for(int j=0;j<numberOfInteractions;j++){
			int index1, coefficient;
			inData>>index1;
			inData>>coefficient;
			hillStruct interStruct;
			interStruct.power=coefficient;
			interStruct.speciesLabel=index1;
			structureVector[j]=interStruct;
			interStruct.hillBool=true;
		}
		inStructure.push_back(structureVector);
	}
}
	




	
	