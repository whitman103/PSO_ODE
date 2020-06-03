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

vector<vector<double> > rungeKuttaIteration(vector<double>& species, vector<vector<hillStruct> > hillIds, double deltaT, double currentTime, double stoppingTime, vector<double> printTimes);
void loadHillStruct(vector<vector<hillStruct> >& inStruct, string toOpen);
double hillFunctionTerm(hillStruct& inStruct,double currentPos);

int main(){
	
	vector<vector<hillStruct> > baseData;
	loadHillStruct(baseData,"trueOutStructure");
	
	vector<double> speciesIn={2,5};
	vector<double> speciesReset=speciesIn;
	
	vector<double> printTimes={2,4,8,16,32};
	
	vector<vector<double> > testingData=rungeKuttaIteration(speciesIn, baseData, 0.01, 0, 65, printTimes);
	
	ofstream trueOutData("trueOutData.txt");
	
	for(int i=0;i<(int)testingData.size();i++){
		trueOutData<<printTimes[i]<<" ";
		for(int j=0;j<(int)testingData[i].size();j++){
			trueOutData<<testingData[i][j]<<" ";
		}
		trueOutData<<endl;
	}
	
	trueOutData.close();
	
	
	speciesIn=speciesReset;
	vector<vector<hillStruct> > testData;
	loadHillStruct(testData,"outFits\\goodFit_26");
	
	testingData=rungeKuttaIteration(speciesIn, testData, 0.01, 0, 65, printTimes);
	
	trueOutData.open("testOutData.txt");
	
	for(int i=0;i<(int)testingData.size();i++){
		trueOutData<<printTimes[i]<<" ";
		for(int j=0;j<(int)testingData[i].size();j++){
			trueOutData<<testingData[i][j]<<" ";
		}
		trueOutData<<endl;
	}
	
	trueOutData.close();
	
	
	
}

void loadHillStruct(vector<vector<hillStruct> >& inStruct, string toOpen){
	ifstream inData(toOpen+".txt");
	if(!inData.good()){
		cout<<"Data not loaded from "+toOpen+".txt"<<endl;
	}
	else{
		string firstIn;
		getline(inData,firstIn);
		int numOfSpecies(0);
		inData>>numOfSpecies;
		cout<<numOfSpecies<<endl;
		inStruct.resize(numOfSpecies);
		int intHold(0);
		double doubleHold(0);
		
		for(int i=0;i<numOfSpecies;i++){
			inData>>intHold; //throwaway label index
			int numOfInteractions(0);
			inData>>numOfInteractions;
			inStruct[i].resize(numOfInteractions);
			for(int j=0;j<numOfInteractions;j++){
				inData>>doubleHold;
				inStruct[i][j].speciesLabel=doubleHold;
				inData>>doubleHold;
				inStruct[i][j].power=doubleHold;
				inData>>doubleHold;
				inStruct[i][j].constant=doubleHold;
				inData>>doubleHold;
				inStruct[i][j].normalization=doubleHold;
				inStruct[i][j].hillBool=true;
			}
		}
	}
	inData.close();
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

double hillFunctionTerm(hillStruct& inStruct,double currentPos){
	if(inStruct.hillBool==true){
		return Hill(inStruct.constant, inStruct.power,currentPos, inStruct.normalization);
	}
	else{
		return invHill(inStruct.constant, inStruct.power, currentPos, inStruct.normalization);
	}
}