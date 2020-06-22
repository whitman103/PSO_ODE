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
#include <filesystem>
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
	
	
	string masterFolder="MasterFolder";
	filesystem::create_directory(masterFolder);
	masterFolder+="//";

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
	
	int numOfRuns(2500);
	double constBound(50), normBound(.05), decayBound(.01);
	vector<double> reportTimes={2,4,8,16};
	//Create "experimental data"
	
	Particle testParticle=Particle(hillDetails,make_tuple(constBound,normBound,decayBound));
	for(int species=0;species<(int)testParticle.normCurrentPos.size();species++){
		for(int interaction=0;interaction<(int)testParticle.normCurrentPos[species].size();interaction++){
			testParticle.normCurrentPos[species][interaction]=randPull()*normBound;
			testParticle.constCurrentPos[species][interaction]=randPull()*constBound;
			testParticle.sampleSolution[species][interaction].normalization=testParticle.normCurrentPos[species][interaction];
			testParticle.sampleSolution[species][interaction].constant=testParticle.constCurrentPos[species][interaction];
		}
	}
	for(int species=0;species<(int)specNum.size();species++){
		testParticle.decayConsts[species]=randPull()*decayBound;
	}
	testParticle.dumpParticleDetails(masterFolder,"testValues");
	
	
	
	vector<vector<vector<double> > > testDist(reportTimes.size(), vector<vector<double> > (specNum.size(),vector<double> (numOfRuns,0)));
	int lastReaction(0);

	for(int run=0;run<numOfRuns;run++){
		specNum=resetSpecNum;
		double runTime(0);
		
		int reactionIndex(0);
		
		for(int i=0;i<(int)testParticle.sampleSolution.size();i++){
			for(int j=0;j<(int)testParticle.sampleSolution[i].size();j++){
				ReactionObject1.reactConsts[reactionIndex]=Hill(testParticle.constCurrentPos[i][j],testParticle.sampleSolution[i][j].power, specNum[testParticle.sampleSolution[i][j].speciesLabel],testParticle.normCurrentPos[i][j]);
				reactionIndex++;
			}
		}
		
		
		
		for(int i=0;i<(int)testParticle.decayConsts.size();i++){
			ReactionObject1.reactConsts[reactionIndex]=testParticle.decayConsts[i];
			reactionIndex++;
		}
		
		int reportIndex(0);
		do{
			reactionIndex=0;
			for(int i=0;i<(int)testParticle.sampleSolution.size();i++){
				for(int j=0;j<(int)testParticle.sampleSolution[i].size();j++){
					if(ReactionObject1.changeCoeffs[lastReaction][testParticle.sampleSolution[i][j].speciesLabel]!=0&&runTime!=0){
						ReactionObject1.reactConsts[reactionIndex]=Hill(testParticle.constCurrentPos[i][j],testParticle.sampleSolution[i][j].power, specNum[testParticle.sampleSolution[i][j].speciesLabel],testParticle.normCurrentPos[i][j]);
					}
					reactionIndex++;
				}
			}
			
			
			tuple<int,double> hold=ReactionObject1.PerformTimeStep2(specNum);
			lastReaction=get<0>(hold);
			runTime+=get<1>(hold);
			ReactionObject1.specChange(specNum,get<0>(hold),ReactionObject1.changeCoeffs);
			
			
			if(runTime>reportTimes[reportIndex]){
				for(int i=0;i<(int)specNum.size();i++){
					testDist[reportIndex][i][run]=specNum[i];
				}
				reportIndex++;
			}
		}while(reportIndex<(int)reportTimes.size());
	}
	
	vector<vector<double> > testMeans(reportTimes.size(), vector<double> (specNum.size(),0));
	for(int i=0;i<(int)reportTimes.size();i++){
		for(int j=0;j<(int)specNum.size();j++){
			testMeans[i][j]=returnMean(testDist[i][j]);
		}
	}
	
	
	
	
	int taskID(-1);
	int nTasks(-1);
	

	
	int numParticles(20);
	vector<Particle> particleSwarm;
	vector<Particle> resetSwarm;
	
	int numInteractionConsts(0);
	
	for(int i=0;i<numParticles;i++){
		Particle interParticle=Particle(hillDetails,make_tuple(constBound,normBound,decayBound));
		for(int species=0;species<(int)interParticle.normCurrentPos.size();species++){
			for(int interaction=0;interaction<(int)interParticle.normCurrentPos[species].size();interaction++){
				interParticle.normCurrentPos[species][interaction]=randPull()*normBound;
				interParticle.constCurrentPos[species][interaction]=randPull()*constBound;
				interParticle.constVelocity[species][interaction]=0;
				interParticle.normVelocity[species][interaction]=0;
				numInteractionConsts+=2;
			}
		}
		interParticle.normBestPos=interParticle.normCurrentPos;
		interParticle.constBestPos=interParticle.constCurrentPos;
		for(int species=0;species<(int)specNum.size();species++){
			interParticle.decayConsts[species]=randPull()*decayBound;
			interParticle.decayVelocities[species]=0;
		}
		interParticle.bestDecayConsts=interParticle.decayConsts;
		particleSwarm.push_back(interParticle);
	}
	numInteractionConsts/=numParticles;
	resetSwarm=particleSwarm;
	
	//Define the search space of the particles
	//Constrain the normalizations to be between -15 and 15
	double inDelta(0);
	/*for(int i=0;i<(int)ReactionObject1.reactConsts.size();i++){
		inDelta+=pow(normBound,2);
	}
	//Constrain the hill constants to be relatively close to the experimental data. 
	for(int i=0;i<(int)ReactionObject1.reactConsts.size();i++){
		inDelta+=pow(constBound,2);
	}
	//Constrain the decay constants to be between infinite and 2 min^{-1}
	for(int i=0;i<(int)specNum.size();i++){
		inDelta+=pow(decayBound,2);
	}*/
	//normalize each parameter direction
	
	
	
	
	double fitnessNormalization(0);
	double bestFitnessValue(1e13);
	double fitnessContainer[numParticles];
	double fitnessValue(0);
	int paramsPerReaction(2);
	
	int sizeOfParameterVector(numInteractionConsts+specNum.size());
	double parameterMatrixHold[(sizeOfParameterVector)*(numParticles)];
	double parameterVectorToSend[sizeOfParameterVector];
	double globalBestParameterSet[sizeOfParameterVector];
	
	inDelta=sqrt(sizeOfParameterVector);
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	FuzzyTree fuzzyStruct(inDelta);
	
	if(nTasks!=numParticles){
		cout<<"Reset Particle Size"<<endl;
		return 0;
	}
	generator.seed(time(NULL)+taskID);
	
	int numOfSets(10);
	for(int set=0;set<numOfSets;set++){
		string currentFolder=masterFolder+"outputSet_"+to_string(set);
		currentFolder+="\\";
		filesystem::create_directory(currentFolder);

		ofstream fitnessValuesOut(currentFolder+"fitnessMasterOut_.txt");
		//reset particles in new positions
		particleSwarm=resetSwarm;

		//first Time setup
		
		Particle myParticle=particleSwarm[taskID];
		
		vector<vector<vector<double> > > outDist(reportTimes.size(), vector<vector<double> > (specNum.size(),vector<double> (numOfRuns,0)));
		lastReaction=0;

		for(int run=0;run<numOfRuns;run++){
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
						if(ReactionObject1.changeCoeffs[lastReaction][myParticle.sampleSolution[i][j].speciesLabel]!=0&&runTime!=0){
							ReactionObject1.reactConsts[reactionIndex]=Hill(myParticle.constCurrentPos[i][j],myParticle.sampleSolution[i][j].power, specNum[myParticle.sampleSolution[i][j].speciesLabel],myParticle.normCurrentPos[i][j]);
						}
						reactionIndex++;
					}
				}
				
				
				tuple<int,double> hold=ReactionObject1.PerformTimeStep2(specNum);
				lastReaction=get<0>(hold);
				runTime+=get<1>(hold);
				ReactionObject1.specChange(specNum,get<0>(hold),ReactionObject1.changeCoeffs);
				
				
				if(runTime>reportTimes[reportIndex]){
					for(int i=0;i<(int)specNum.size();i++){
						outDist[reportIndex][i][run]=specNum[i];
					}
					reportIndex++;
				}
			}while(reportIndex<(int)reportTimes.size());
		}
		
		vector<vector<double> > outMeans(reportTimes.size(), vector<double> (specNum.size(),0));
		for(int i=0;i<(int)reportTimes.size();i++){
			for(int j=0;j<(int)specNum.size();j++){
				outMeans[i][j]=returnMean(outDist[i][j]);
			}
		}
		myParticle.lastFitness=0;
		fitnessValue=0;
		for(int i=0;i<(int)testMeans.size();i++){
			for(int j=0;j<(int)testMeans[i].size();j++){
				fitnessValue+=pow(outMeans[i][j]-testMeans[i][j],2);
			}
		}
		myParticle.currentFitness=fitnessValue;
		fitnessNormalization=fitnessValue;
		
		if(myParticle.currentFitness<myParticle.bestFitness){
			myParticle.normBestPos=myParticle.normCurrentPos;
			myParticle.constBestPos=myParticle.constCurrentPos;
			myParticle.bestDecayConsts=myParticle.decayConsts;
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

		
		MPI_Gather(&myParticle.currentFitness, 1, MPI_DOUBLE, fitnessContainer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		if(taskID==0){
			int bestParticle(0);
			double worstFitness(0);
			bool setNewGlobalBest(false);
			for(int i=0;i<numParticles;i++){
				if(fitnessContainer[i]<bestFitnessValue){
					bestFitnessValue=fitnessContainer[i];
					setNewGlobalBest=true;
					bestParticle=i;
				}
				if(fitnessContainer[i]>worstFitness){
					worstFitness=fitnessContainer[i];
				}
			}
			if(setNewGlobalBest){
				for(int i=0;i<sizeOfParameterVector;i++){
					globalBestParameterSet[i]=parameterMatrixHold[bestParticle*sizeOfParameterVector+i];
				}
			}
			fitnessNormalization=worstFitness;
			for(int i=0;i<sizeOfParameterVector;i++){
				parameterVectorToSend[i]=globalBestParameterSet[i];
			}
		}
		MPI_Bcast(&fitnessNormalization, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		fuzzyStruct.phiNormalization=fitnessNormalization;
		MPI_Bcast(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//first update occurs with mean hyperparameters by fiat
		fuzzyStruct.delta=sqrt(myParticle.performUpdate(&generator, &fuzzyStruct, parameterVectorToSend));
		
		
		
		// Iterate the solutions for PSO
		int numIterations(100);
		for(int iteration=0;iteration<numIterations;iteration++){
			//Generate distributions
			for(int run=0;run<numOfRuns;run++){
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
							if(ReactionObject1.changeCoeffs[lastReaction][myParticle.sampleSolution[i][j].speciesLabel]!=0&&runTime!=0){
								ReactionObject1.reactConsts[reactionIndex]=Hill(myParticle.constCurrentPos[i][j],myParticle.sampleSolution[i][j].power, specNum[myParticle.sampleSolution[i][j].speciesLabel],myParticle.normCurrentPos[i][j]);
							}
							reactionIndex++;
						}
					}
					
					
					tuple<int,double> hold=ReactionObject1.PerformTimeStep2(specNum);
					lastReaction=get<0>(hold);
					runTime+=get<1>(hold);
					ReactionObject1.specChange(specNum,get<0>(hold),ReactionObject1.changeCoeffs);
					
					
					if(runTime>reportTimes[reportIndex]){
						for(int i=0;i<(int)specNum.size();i++){
							outDist[reportIndex][i][run]=specNum[i];
						}
						reportIndex++;
					}
				}while(reportIndex<(int)reportTimes.size());
			}
			
			for(int i=0;i<(int)reportTimes.size();i++){
				for(int j=0;j<(int)specNum.size();j++){
					outMeans[i][j]=returnMean(outDist[i][j]);
				}
			}
			myParticle.lastFitness=myParticle.currentFitness;
			fitnessValue=0;
			for(int i=0;i<(int)testMeans.size();i++){
				for(int j=0;j<(int)testMeans[i].size();j++){
					fitnessValue+=pow(outMeans[i][j]-testMeans[i][j],2);
				}
			}
			myParticle.currentFitness=fitnessValue;
			if(myParticle.currentFitness<myParticle.bestFitness){
				myParticle.normBestPos=myParticle.normCurrentPos;
				myParticle.constBestPos=myParticle.constCurrentPos;
				myParticle.bestDecayConsts=myParticle.decayConsts;
				myParticle.bestFitness=myParticle.currentFitness;
			}
			fuzzyStruct.setParameters();
			
			fillingIndex=0;
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
			
			MPI_Gather(&myParticle.currentFitness, 1, MPI_DOUBLE, fitnessContainer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			int bestParticle(0);
			if(taskID==0){
				bool setNewGlobalBest(false);
				for(int i=0;i<numParticles;i++){
					fitnessValuesOut<<fitnessContainer[i]<<" ";
					if(fitnessContainer[i]<bestFitnessValue){
						bestFitnessValue=fitnessContainer[i];
						setNewGlobalBest=true;
						bestParticle=i;
					}
					if(fitnessContainer[i]>fitnessNormalization){
						fitnessNormalization=fitnessContainer[i];
					}
				}
				fitnessValuesOut<<endl;
				
				if(setNewGlobalBest){
					for(int i=0;i<sizeOfParameterVector;i++){
						globalBestParameterSet[i]=parameterMatrixHold[bestParticle*sizeOfParameterVector+i];
					}
				}
				for(int i=0;i<sizeOfParameterVector;i++){
					parameterVectorToSend[i]=globalBestParameterSet[i];
				}
			}
			MPI_Bcast(&bestParticle, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&fitnessNormalization, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			fuzzyStruct.phiNormalization=fitnessNormalization;
			MPI_Bcast(parameterVectorToSend, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			fuzzyStruct.calculatePhi(myParticle.lastFitness,myParticle.currentFitness);
			fuzzyStruct.delta=sqrt(myParticle.performUpdate(&generator, &fuzzyStruct, parameterVectorToSend));
		}
		
		

		for(int i=0;i<(int)myParticle.sampleSolution.size();i++){
			for(int j=0;j<(int)myParticle.sampleSolution[i].size();j++){
				myParticle.sampleSolution[i][j].constant=myParticle.constBestPos[i][j];
				myParticle.sampleSolution[i][j].normalization=myParticle.normBestPos[i][j];
			}
		}
		
		
		myParticle.dumpParticleDetails(currentFolder,to_string(taskID));
		if(taskID==0){
			fillingIndex=0;
			for(int i=0;i<(int)myParticle.sampleSolution.size();i++){
				for(int j=0;j<(int)myParticle.sampleSolution[i].size();j++){
					myParticle.sampleSolution[i][j].normalization=globalBestParameterSet[fillingIndex];
					fillingIndex++;
				}
			}
			for(int i=0;i<(int)myParticle.sampleSolution.size();i++){
				for(int j=0;j<(int)myParticle.sampleSolution[i].size();j++){
					myParticle.sampleSolution[i][j].constant=globalBestParameterSet[fillingIndex];
					fillingIndex++;
				}
			}
			for(int i=0;i<(int)myParticle.decayConsts.size();i++){
				myParticle.decayConsts[i]=globalBestParameterSet[fillingIndex];
				fillingIndex++;
			}
			myParticle.dumpParticleDetails(currentFolder,"globalBest");
		}
		fitnessValuesOut.close();
	}
			
	
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
			if(outStruct[i][j].hillBool==true){
				StructCoeffs<<1<<" "<<i<<" "<<-1<<endl;
			}
			else{
				StructCoeffs<<1<<" "<<i<<" "<<1<<endl;
			}
		}
	}
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<1<<" "<<i<<" "<<1<<endl;
	}
	
	StructCoeffs.close();
	StructCoeffs.open(outLeader+"Consts.txt");
	
	StructCoeffs<<outStruct.size()<<endl;
	for(int i=0;i<(int)outStruct.size();i++){
		StructCoeffs<<randSite(125)<<endl;
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
			interStruct.hillBool=randPull()>0.25;
			structureVector[j]=interStruct;
		}
		inStructure.push_back(structureVector);
	}
}
	




	
	