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
#include "fuzzyDef.h"
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
	vector<double> reportTimes={2,4,8,16,32,64};
	
	int numParticles(2);
	vector<Particle> particleSwarm;
	
	double constBound(0), normBound(0), decayBound(0);
	
	for(int i=0;i<numParticles;i++){
		Particle interParticle=Particle(hillDetails,make_tuple(constBound,normBound,decayBound));
		for(int species=0;species<(int)interParticle.normCurrentPos.size();species++){
			for(int interaction=0;interaction<(int)interParticle.normCurrentPos[species].size();interaction++){
				interParticle.normCurrentPos[species][interaction]=twoSidedRandPull()*normBound;
				interParticle.constCurrentPos[species][interaction]=randPull()*constBound;
				
			}
		}
		interParticle.normBestPos=interParticle.normCurrentPos;
		interParticle.constBestPos=interParticle.constCurrentPos;
		for(int species=0;species<(int)specNum.size();species++){
			interParticle.decayConsts[species]=randPull()*decayBound;
		}
		interParticle.bestDecayConsts=interParticle.decayConsts;
		particleSwarm.push_back(interParticle);
	}
	
	//Define the search space of the particles
	//Constrain the normalizations to be between -15 and 15
	double inDelta(0);
	for(int i=0;i<(int)ReactionObject1.reactConsts.size();i++){
		inDelta+=pow(2*normBound,2);
	}
	//Constrain the hill constants to be relatively close to the experimental data. 
	for(int i=0;i<(int)ReactionObject1.reactConsts.size();i++){
		inDelta+=pow(constBound,2);
	}
	//Constrain the decay constants to be between infinite and 2 min^{-1}
	for(int i=0;i<(int)specNum.size();i++){
		inDelta+=pow(decayBound,2);
	}
	inDelta=sqrt(inDelta);
	
	
	
	
	double bestFitnessValue(0);
	double fitnessContainer[numOfParticles];
	double fitnessValue(0);
	int paramsPerReaction(2);
	
	int sizeOfParameterVector(paramsPerReaction*ReactionObject1.reactConsts.size()+specNum.size());
	double parameterMatrixHold[(paramsPerReaction*ReactionObject1.reactConsts.size()+specNum.size())*(numParticles+1)];
	double parameterVectorToSend[(paramsPerReaction*ReactionObject1.reactConsts.size()+specNum.size())];
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	FuzzyTree fuzzyStruct(inDelta);
	
	if(nTasks!=numParticles){
		cout<<"Reset Particle Size"<<endl;
		return 0;
	}
	generator.seed(time(NULL)+taskID);
	
	
	
			
	//first Time setup
	
	Particle myParticle=particleSwarm[taskID];
	
	vector<vector<vector<double> > > outDist(reportTimes, vector<vector<double> > (specNum.size(),vector<double> (numOfRuns,0)));
	int lastReactionChange(0);

	for(int run=0;i<numOfRuns;run++){
		specNum=resetSpecNum;
		double runTime(0);
		
		int reactionIndex(0);
		
		for(int i=0;i<(int)myParticle.sampleSolution.size();i++){
			for(int j=0;j<(int)myParticle.sampleSolution[i].size();j++){
				ReactionObject1.reactConsts[reactionIndex]=Hill(myParticle.constCurrentPos[i][j],myParticle.sampleSolution[i][j].power, specNum[myParticle.sampleSolution[i][j].speciesLabel],myParticle.normCurrentPos[i][j]);
					reactionIndex++;
			}
		}
		
		for(int i=0;i<(int)myParticle.decayConsts.size();i++){
			ReactionObject1.reactConsts[reactionIndex]=myParticle.decayConsts[i];
			reactionIndex++;
		}
		
		int reportIndex(0);
		do{
			reactionIndex=0;
			
			for(int i=0;i<(int)myParticle.sampleSolution.size();i++){
				for(int j=0;j<(int)myParticle.sampleSolution[i].size();j++){
					if(myParticle.sampleSolutions[i][j].speciesLabel==lastReactionChange&&runTime!=0){
						ReactionObject1.reactConsts[reactionIndex]=Hill(myParticle.constCurrentPos[i][j],myParticle.sampleSolution[i][j].power, specNum[myParticle.sampleSolution[i][j].speciesLabel],myParticle.normCurrentPos[i][j]);
					}
					reactionIndex++;
				}
			}
				
			for(int i=0;i<(int)myParticle.decayConsts.size();i++){
				ReactionObject1.reactConsts[reactionIndex]=myParticle.decayConsts[i];
				reactionIndex++;
			}
			
			
			tuple<int,double> hold=ReactionObject1.PerformTimeStep2(specNum);
			runTime+=get<1>(hold);
			ReactionObject1.specChange(specNum,get<0>(hold),ReactionObject1.changeCoeffs);
			for(int i=0;i<(int)ReactionObject1[get<0>(hold)].changeCoeffs.size();i++){
				if(ReactionObject1[get<0>(hold)][i]!=0){
					lastReactionChange=i;
				}
			}
			if(runTime>reportTimes[reportIndex]){
				for(int i=0;i<(int)specNum.size();i++){
					outDist[reportIndex][i][run]=specNum[i];
				}
				reportIndex++;
			}
		}while(reportIndex<reportTimes.size());
	}
	
	vector<vector<double> > testMeans(reportTimes.size(), vector<double> (numSpecies.size(),0));
	for(int i=0;i<(int)reportTimes.size();i++){
		for(int j=0;j<(int)numSpecies.size();j++){
			testMeans[i][j]=returnMean(outDist[i][j]);
		}
	}
	
	double fitnessValue(0);
	for(int i=0;i<(int)testMeans.size();i++){
		for(int j=0;j<(int)testMeans[i].size();j++){
			fitnessValue+=pow(testMeans[i][j]-experimentalMeans[i][j],2);
		}
	}
	
	int fillingIndex(0);
	for(int i=0;i<(int)myParticle.normCurrentPos.size();i++){
		for(int j=0;j<(int)myParticle.normCurrentPos[i].size();j++){
			parameterVectorToSend[fillingIndex]=myParticle.normCurrentPos[i][j];
			fillingIndex++;
		}
	}
	for(int i=0;i<(int)myParticle.constCurrentPos.size();i++){
		for(int j=0;j<(int)myParticle.constCurrentPos[i].size();j++){
			parameterVectorToSend[fillingIndex]=myParticle.constCurrentPos[i][j];
			fillingIndex++;
		}
	}
	for(int i=0;i<(int)myParticle.decayConsts.size();i++){
		parameterVectorToSend[fillingIndex]=myParticle.decayConsts[i];
		fillingIndex++;
	}

	
	MPI_Gather(&fitnessValue, 1, MPI_DOUBLE, fitnessContainer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(taskID==0){
		int bestParticle(0);
		for(int i=0;i<numParticles;i++){
			if(fitnessContainer[i]<masterFitnessValue){
				bestFitnessValue=fitnessContainer[i];
				bestParticle=i;
			}
		}
		for(int i=0;i<sizeOfParameterVector;i++){
			parameterVectorToSend[i]=parameterMatrixHold[bestParticle*sizeOfParameterVector+i];
		}
	}
	MPI_Bcast(&bestFitnessValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

	//first update occurs with mean hyperparameters by fiat
	int fillingIndex(0);
	fuzzyStruct.delta=0;
	for(int i=0;i<(int)myParticle.normCurrentPos.size();i++){
		for(int j=0;j<(int)myParticle.normCurrentPos[i].size();j++){
			double rand1(randPull());
			double rand2(randPull());
			double proposedVelocity(0);
			proposedVelocity+=fuzzyStruct.inertia*myParticle.normCurrentPos[i][j];
			proposedVelocity+=fuzzyStruct.social*rand1*(myParticle.normCurrentPos[i][j]-parameterVectorToSend[fillingIndex]);
			proposedVelocity+=fuzzyStruct.cognitive*rand2*(myParticle.normBestPos[i][j]);
			fuzzyStruct.delta+=pow(proposedVelocity,2);
			//Checks velocity limits
			if(proposedUpdate<fuzzyStruct.U*(-2.*normBound)){
				proposedUpdate=fuzzyStruct.U*(-2.*normBound);
			}
			if(proposedUpdate>fuzzyStruct.U*(2.*normBound)){
				proposedUpdate=fuzzyStruct.U*(2.*normBound);
			}
			//Checks bounds of parameter space;
			if(myParticle.normCurrentPos[i][j]+proposedUpdate<(-1.*normBound)){
				myParticle.normCurrentPos[i][j]=-1.*normBound+(-1.*proposedUpdate);
			}
			else{
				if(myParticle.normCurrentPos[i][j]+proposedUpdate>normBound){
					myParticle.normCurrentPos[i][j]=normBound-proposedUpdate;
				}
				else{
					myParticle.normCurrentPos[i][j]+=proposedUpdate;
				}
			}
			fillingIndex++;
		}
	}
	
	
	for(int i=0;i<(int)myParticle.constCurrentPos.size();i++){
		for(int j=0;j<(int)myParticle.normCurrentPos[i].size();j++){
			double rand1(randPull());
			double rand2(randPull());
			double proposedVelocity(0);
			proposedVelocity+=fuzzyStruct.inertia*myParticle.constCurrentPos[i][j];
			proposedVelocity+=fuzzyStruct.social*rand1*(myParticle.constCurrentPos[i][j]-parameterVectorToSend[fillingIndex]);
			proposedVelocity+=fuzzyStruct.cognitive*rand2*(myParticle.normBestPos[i][j]);
			fuzzyStruct.delta+=pow(proposedVelocity,2);
			//Checks velocity limits
			if(proposedUpdate<fuzzyStruct.U*(-1*.constBound){
				proposedUpdate=fuzzyStruct.U*(constBound);
			}
			if(proposedUpdate>fuzzyStruct.U*(constBound){
				proposedUpdate=fuzzyStruct.U*(constBound);
			}
			//Checks bounds of parameter space;
			if(myParticle.constCurrentPos[i][j]+proposedUpdate<0){
				myParticle.constCurrentPos[i][j]=(-1.*proposedUpdate);
			}
			else{
				if(myParticle.constCurrentPos[i][j]+proposedUpdate>normBound){
					myParticle.constCurrentPos[i][j]=normBound-proposedUpdate;
				}
				else{
					myParticle.constCurrentPos[i][j]+=proposedUpdate;
				}
			}
			fillingIndex++;
		}
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
		outHold+=inList[i];
	}
	outHold/=(double)inList.size();
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
	




	
	