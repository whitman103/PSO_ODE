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

inline double invHill(double constant, double power, double argument, double norm){
	return norm*pow(constant,power)/(pow(constant,power)+pow(argument,power));
}

inline double Hill(double constant, double power, double argument,double norm){
	return norm*pow(argument,power)/(pow(constant,power)+pow(argument,power));
}


/* Custom structures */

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
	double bestWellness;
	double currentWellness;
} particle;


/* Function declarations */
	
void loadExperimentData(vector<double>& inData);
double hillFunctionTerm(hillStruct& inStruct,double currentPos);
vector<vector<double> > rungeKuttaIteration(vector<double>& species, vector<vector<hillStruct> > hillIds, double deltaT, double currentTime, double stoppingTime, vector<double> printTimes);
void loadHillStructDetails(vector<vector<hillStruct> >& inStructure, string inFile);
void dumpStructure(particle inParticle, int particleCount);
	
int main(){
	
	generator.seed(time(NULL));
	//system("./cpExe.sh");
	//system("./spk_serial < in.Base");
	
	
	
	
	
	string baseHillStructure("randomHillStructure");
	
	vector<vector<hillStruct> > hillDetails;
	loadHillStructDetails(hillDetails, baseHillStructure);
	
	int speciesInSize(hillDetails.size());
	vector<double> speciesIn={2,5,10,3,5};
	vector<double> speciesReset=speciesIn;
	
	vector<vector<tuple<double,double> > > baseParameters(speciesIn.size());
	
	for(int i=0;i<(int)hillDetails.size();i++){
		for(int j=0;j<(int)hillDetails[i].size();j++){
			hillDetails[i][j].hillBool=true;
			hillDetails[i][j].normalization=3.*randPull();
			hillDetails[i][j].constant=randSite(15);
			baseParameters[i].push_back(make_tuple(hillDetails[i][j].normalization,hillDetails[i][j].constant));
		}
	}
	
	ofstream trueOutStructure("trueOutStructure.txt");
	
	trueOutStructure<<"#m N: species Power Constant Norm"<<endl;
	trueOutStructure<<hillDetails.size()<<endl;
	for(int i=0;i<(int)hillDetails.size();i++){
		trueOutStructure<<i<<" "<<hillDetails[i].size()<<" ";
		for(int j=0;j<(int)hillDetails[i].size();j++){
			trueOutStructure<<hillDetails[i][j].speciesLabel<<" "<<hillDetails[i][j].power<<" "<<hillDetails[i][j].constant<<" "<<hillDetails[i][j].normalization<<" ";
		}
		trueOutStructure<<endl;
	}
	trueOutStructure.close();
	
	vector<double> printTimes={2,4,8,16,32};
	
	vector<vector<double> > testingData=rungeKuttaIteration(speciesIn, hillDetails, 0.01, 0, 65, printTimes);
	
	ofstream trueOutData("trueOutData.txt");
	
	for(int i=0;i<(int)testingData.size();i++){
		trueOutData<<printTimes[i]<<" ";
		for(int j=0;j<(int)testingData[i].size();j++){
			trueOutData<<testingData[i][j]<<" ";
		}
		trueOutData<<endl;
	}
	
	trueOutData.close();
	
	const int numParticles(20);
	
	vector<particle> particleSwarm(numParticles);
	
	//initialize particleSwarm
	
	for(int i=0;i<(int)particleSwarm.size();i++){
		particle interParticle;
		loadHillStructDetails(interParticle.sampleSolution,baseHillStructure);
		interParticle.currentPosition.resize(interParticle.sampleSolution.size());
		interParticle.currentVelocity.resize(interParticle.sampleSolution.size());
		interParticle.currentWellness=0;
		for(int species=0;species<(int)interParticle.sampleSolution.size();species++){
			interParticle.currentPosition[species].resize(interParticle.sampleSolution[species].size());
			interParticle.currentVelocity[species].resize(interParticle.sampleSolution[species].size());
			for(int interaction=0;interaction<(int)interParticle.sampleSolution[species].size();interaction++){
				interParticle.currentPosition[species][interaction]=make_tuple(15*randPull(),randSite(30));
				interParticle.currentVelocity[species][interaction]=make_tuple(0,0);
			}
		}
		interParticle.bestPosition=interParticle.currentPosition;
		particleSwarm[i]=interParticle;
	}
	
	// calculate first set of fits
	
	int indexOfBestFit(0);
	double currentBestMSE(1e13);
	for(int sol=0;sol<(int)particleSwarm.size();sol++){
		particle* currentParticle=&particleSwarm[sol];
		for(int i=0;i<(int)(*currentParticle).sampleSolution.size();i++){
			for(int j=0;j<(int)(*currentParticle).sampleSolution[i].size();j++){
				auto[interNorm,interConstant]=(*currentParticle).currentPosition[i][j];
				(*currentParticle).sampleSolution[i][j].normalization=interNorm;
				(*currentParticle).sampleSolution[i][j].constant=interConstant;
				(*currentParticle).sampleSolution[i][j].hillBool=true;
			}
		}
		speciesIn=speciesReset;
		vector<vector<double> > currentTest=rungeKuttaIteration(speciesIn, (*currentParticle).sampleSolution, 0.01, 0, 65, printTimes);
			
		double MSE(0);
		for(int i=0;i<(int)currentTest.size();i++){
			for(int j=0;j<(int)currentTest[i].size();j++){
				MSE+=pow(currentTest[i][j]-testingData[i][j],2);
			}
		}
		(*currentParticle).currentWellness=MSE;
		(*currentParticle).bestWellness=MSE;
		if(MSE<currentBestMSE){
			indexOfBestFit=sol;
			currentBestMSE=MSE;
		}
	}
	
	cout<<currentBestMSE<<endl;
	//first update of solutions
	int goodSolutions(0);
	particle bestParticle;
	for(int iteration=0;iteration<1000;iteration++){
		cout<<iteration<<endl;
		cout<<currentBestMSE<<endl;
		for(int sol=0;sol<(int)particleSwarm.size();sol++){
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
		}
	
	
	
		
		for(int sol=0;sol<(int)particleSwarm.size();sol++){
			particle* currentParticle=&particleSwarm[sol];
			for(int i=0;i<(int)(*currentParticle).sampleSolution.size();i++){
				for(int j=0;j<(int)(*currentParticle).sampleSolution[i].size();j++){
					auto[interNorm,interConstant]=(*currentParticle).currentPosition[i][j];
					(*currentParticle).sampleSolution[i][j].normalization=interNorm;
					(*currentParticle).sampleSolution[i][j].constant=interConstant;
				}
			}
			speciesIn=speciesReset;
			vector<vector<double> > currentTest=rungeKuttaIteration(speciesIn, (*currentParticle).sampleSolution, 0.01, 0, 65, printTimes);
			
			
			double MSE(0);
			for(int i=0;i<(int)currentTest.size();i++){
				for(int j=0;j<(int)currentTest[i].size();j++){
					MSE+=pow(currentTest[i][j]-testingData[i][j],2);
				}
			}
			(*currentParticle).currentWellness=MSE;
			if(MSE<currentBestMSE){
				indexOfBestFit=sol;
				currentBestMSE=MSE;
				bestParticle=(*currentParticle);
			}
			if(MSE<(*currentParticle).bestWellness){
				(*currentParticle).bestPosition=(*currentParticle).currentPosition;
				(*currentParticle).bestWellness=MSE;
			}
			if(MSE<1000){
				dumpStructure((*currentParticle),goodSolutions);
				goodSolutions++;
			}
		}
		cout<<indexOfBestFit<<" "<<currentBestMSE<<endl;
		cout<<endl;
		
	}
	
	for(int i=0;i<(int)particleSwarm.size();i++){
		dumpStructure(particleSwarm[i], goodSolutions);
		goodSolutions++;
	}
	
	dumpStructure(bestParticle,goodSolutions);
	
	cout<<currentBestMSE<<endl;
	
	
	
	return 0;
}

void dumpStructure(particle inParticle, int particleCount){
	ofstream outParticle("outFits\\goodFit_"+to_string(particleCount)+".txt");
	if(!outParticle.good()){
		cout<<"hold up"<<endl;
	}
	outParticle<<"#m N: species Power Constant Norm"<<endl;
	outParticle<<inParticle.sampleSolution.size()<<endl;
	
	for(int i=0;i<(int)inParticle.sampleSolution.size();i++){
		outParticle<<i<<" "<<inParticle.sampleSolution[i].size()<<" ";
		for(int j=0;j<(int)inParticle.sampleSolution[i].size();j++){
			outParticle<<inParticle.sampleSolution[i][j].speciesLabel<<" "<<inParticle.sampleSolution[i][j].power<<" "<<inParticle.sampleSolution[i][j].constant<<" "<<inParticle.sampleSolution[i][j].normalization<<" ";
		}
		outParticle<<endl;
	}
	outParticle.close();
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

/* Load in synthetic or true data distributions for comparison and loss function calculations. Returns a normalized distribution in inDist. */
void loadExperimentData(vector<double>& inDist){
	vector<int> inData;
	ifstream myfile("outFile.txt");
	int inSize(0);
	myfile>>inSize;
	for(int i=0;i<inSize;i++){
		int dataHold(0);
		for(int j=0;j<4;j++){
			myfile>>dataHold;
		}
		inData.push_back(dataHold);
		myfile>>dataHold;
	}
	myfile.close();
	
	int distSize=*max_element(inData.begin(),inData.end());
	inDist.resize(distSize+1,0);
	
	for(int i=0;i<(int)inData.size();i++){
		inDist[inData[i]]+=1./(double)inData.size();
	}
}

/* Quick function for calculating Hill functions using the structure defined in the preamble. Constituent function calls are enumerated by the compiler for speed. */

double hillFunctionTerm(hillStruct& inStruct,double currentPos){
	if(inStruct.hillBool==true){
		return Hill(inStruct.constant, inStruct.power,currentPos, inStruct.normalization);
	}
	else{
		return invHill(inStruct.constant, inStruct.power, currentPos, inStruct.normalization);
	}
}

vector<vector<double> > rungeKuttaIteration(vector<double>& species, vector<vector<hillStruct> > hillIds, double deltaT, double currentTime, double stoppingTime, vector<double> printTimes){
	ofstream outSubstances("testOut_species.txt");
	vector<vector<double> > outputSet;
	int n=(int)(((stoppingTime-currentTime))/(deltaT));
	int printIndex(0);
	printTimes.push_back(stoppingTime+100);
	for(int t=0;t<n;t++){
	
		vector<double> k1(species.size(),0);
		for(int i=0;i<(int)k1.size();i++){
			for(int j=0;j<(int)hillIds[i].size();j++){
				k1[i]+=deltaT*hillFunctionTerm(hillIds[i][j],species[hillIds[i][j].speciesLabel]);
			}
		}
		vector<double> k2(species.size(),0);
		for(int i=0;i<(int)k2.size();i++){
			for(int j=0;j<(int)hillIds[i].size();j++){
				k2[i]+=deltaT*hillFunctionTerm(hillIds[i][j],species[hillIds[i][j].speciesLabel]+k1[i]/2.);
			}
		}
		vector<double> k3(species.size(),0);
		for(int i=0;i<(int)k3.size();i++){
			for(int j=0;j<(int)hillIds[i].size();j++){
				k3[i]+=deltaT*hillFunctionTerm(hillIds[i][j],species[hillIds[i][j].speciesLabel]+k2[i]/2.);
			}
		}
		vector<double> k4(species.size(),0);
		for(int i=0;i<(int)k4.size();i++){
			for(int j=0;j<(int)hillIds[i].size();j++){
				k4[i]+=deltaT*hillFunctionTerm(hillIds[i][j],species[hillIds[i][j].speciesLabel]+k3[i]);
			}
		}
		for(int i=0;i<(int)species.size();i++){
			species[i]=species[i]+k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.;
		}
		if(t*deltaT>printTimes[printIndex]){
			outputSet.push_back(species);
			printIndex++;
		}
	}
	
	outSubstances.close();
	return outputSet;
}