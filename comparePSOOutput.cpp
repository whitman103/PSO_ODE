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

using namespace std;

int main(){


    ifstream inData("MasterFolder\\particleDump_testValues.txt");
    string throwData("");
    double inHold(0);
    getline(inData,throwData);
    int numSpecies(0), numPlots(0), hillParameters(2);
   
    inData>>numSpecies;
    inData>>inHold;
    numPlots=inHold;
    //trueParameters[0] are Constants, trueParameters[1] are norms
    vector<vector<double> > trueParameters(numPlots,vector<double>(hillParameters,0));
    vector<double> trueDecay(numSpecies,0);

    int fillingIndex(0);
    for(int i=0;i<numSpecies;i++){
        //Throw out first input, simply labels species
        inData>>inHold;
        int currentNumberReactions(0);
        inData>>currentNumberReactions;
        for(int j=0;j<currentNumberReactions;j++){
            inData>>inHold;//species Label
            inData>>inHold;//power 
            inData>>trueParameters[fillingIndex][0];//Constant
            inData>>trueParameters[fillingIndex][1];//Norm
            inData>>inHold;//Pos/Neg
            fillingIndex++;
        }
    }
    for(int i=0;i<numSpecies;i++){
        inData>>trueDecay[i];
    }
    inData.close();

    

    ofstream outPlot("trueOutPlot.txt");
    for(int i=0;i<(int)trueParameters.size();i++){
        for(int j=0;j<(int)trueParameters[i].size();j++){
            outPlot<<trueParameters[i][j]<<" ";
        }
    }
    for(int i=0;i<numSpecies;i++){
        outPlot<<trueDecay[i]<<" ";
    }
    outPlot.close();

    outPlot.open("fitData.txt");

    for(int fileIndex=0;fileIndex<5;fileIndex++){
        inData.open("MasterFolder\\particleDump_globalBest_"+to_string(fileIndex)+".txt");
        getline(inData,throwData);
        inData>>numSpecies;//Throw extras
        inData>>inHold;//Throw extras
        if(!inData.good()){
            cout<<"particle doesn't exist";
            return 0;
        }
        fillingIndex=0;
        for(int i=0;i<numSpecies;i++){
            //Throw out first input, simply labels species
            inData>>inHold;
            int currentNumberReactions(0);
            inData>>currentNumberReactions;
            for(int j=0;j<currentNumberReactions;j++){
                inData>>inHold;//species Label
                inData>>inHold;//power 
                inData>>trueParameters[fillingIndex][0];//Constant
                inData>>trueParameters[fillingIndex][1];//Norm
                inData>>inHold;//Pos/Neg
                fillingIndex++;
            }
        }
        for(int i=0;i<numSpecies;i++){
            inData>>trueDecay[i];
        }
        for(int i=0;i<(int)trueParameters.size();i++){
            for(int j=0;j<(int)trueParameters[i].size();j++){
                outPlot<<trueParameters[i][j]<<" ";
            }
        }
        for(int i=0;i<numSpecies;i++){
            outPlot<<trueDecay[i]<<" ";
        }
        outPlot<<endl;
        inData.close();
    }

    outPlot.close();


    return 0;
}