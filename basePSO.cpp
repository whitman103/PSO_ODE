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
	hillStruct sampleSolution;
	vector<vector<tuple<double,double> > currentPosition;
	vector<vector<tuple<double,double> > currentVelocity;
	double currentWellness;
} particle;


/* Function declarations */
	
void loadExperimentData(vector<double>& inData);
double hillFunctionTerm(hillStruct& inStruct,double currentPos);
void rungeKuttaIteration(vector<double>& species, vector<vector<hillStruct> > hillIds, double deltaT, double currentTime, double stoppingTime);
void loadHillStructDetails(vector<vector<HillStruct> >& inStructure);
	
int main(){
	
	generator.seed(time(NULL));
	system("./cpExe.sh");
	system("./spk_serial < in.Base");
	
	
	vector<double> experDist;
	loadExperimentData(experDist);
	
	vector<double> currentBestParameters={2,.1,4,0.0005,.6,.2};
	
	
	ofstream nextParameters("candidateParameterSets\\nextParameters.txt");
	int swarmSize(10);
	nextParameters<<swarmSize<<"\n";
	for(int i=0;i<swarmSize;i++){
		vector<double> prospectiveParameters=currentBestParameters;
		prospectiveParameters[randSite(prospectiveParameters.size())]*=(1+.1*(2*randPull()-1.));
		for(int j=0;j<(int)prospectiveParameters.size();j++){
			nextParameters<<prospectiveParameters[j]<<",";
		}
		nextParameters<<endl;
	}
	nextParameters.close();
	system("python createSpparksInputs.py");
	
	
	const int numOfSpecies(2);
	
	vector<double> speciesIn={2,5};
	
	vector<vector<hillStruct> > hillDetails;
	loadHillStructDetails(hillDetails);
	
	vector<vector<tuple<double,double> > > baseParameters(speciesIn.size());
	
	for(int i=0;i<(int)hillDetails.size();i++){
		for(int j=0;j<(int)hillDetails[i].size();j++){
			hillDetails[i][j].hillBool=true;
			hillDetails[i][j].normalization=3.*randPull();
			hillDetails[i][j].constant=randSite(15);
			baseParameters[i].push_back(make_tuple(hillDetails[i][j].normalization,hillDetails[i][j].constant));
		}
	}
	
	vector<double> printTimes={2,4,8,16,32};
	
	vector<vector<double> > testingData=rungeKuttaIteration(speciesIn, hillDetails, 0.01, 0, 100, printTimes);
	
	
	
	
	
	
	return 0;
}


void loadHillStructDetails(vector<vector<HillStruct> >& inStructure){
	ifstream inData("hillStructureDefintion.txt");
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
			interStruct.pow=coefficient;
			interStruct.speciesLabel=index1;
			structureVector[j]=interStruct;
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
	cout<<n<<endl;
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