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

    string inputFolder("ThreeFolder");
    ifstream inData(inputFolder+"//particleDump_testValues.txt");
    if(!inData.good()){
        cout<<"can't load data";
        return 0;
    }
    string throwData("");
    double inHold(0);
    getline(inData,throwData);
    int numSpecies(0), numPlots(0), hillParameters(3), numParticles(25);
   
    ofstream texOut("texOutput.tex");
    texOut<<"\\begin{equation}"<<endl;
    texOut<<"\\begin{aligned}"<<endl;

    inData>>numSpecies;
    inData>>inHold;
    numPlots=inHold;
    //trueParameters[0] are Constants, trueParameters[1] are norms, trueParameters[2] are powers
    vector<vector<double> > trueParameters(numPlots,vector<double>(hillParameters,0));
    vector<double> trueDecay(numSpecies,0);

    int fillingIndex(0);
    for(int i=0;i<numSpecies;i++){
        texOut<<"\\pd{\\langle N_"+to_string(i)+"\\rangle}{t}&=\\left(";
        bool inverseOrNot(false);
        int speciesLabel(0);
        //Throw out first input, simply labels species
        inData>>inHold;
        int currentNumberReactions(0);
        inData>>currentNumberReactions;
        for(int j=0;j<currentNumberReactions;j++){
            inData>>speciesLabel;//species Label
            inData>>trueParameters[fillingIndex][2];//power 
            inData>>trueParameters[fillingIndex][0];//Constant
            inData>>trueParameters[fillingIndex][1];//Norm
            inData>>inHold;//Pos/Neg
            if(inHold==0){
               texOut<<"-d_{"+to_string(fillingIndex)+"}\\frac{N_"+to_string(speciesLabel)+"^2}{N_"+to_string(speciesLabel)+"^2+K_{"+to_string(fillingIndex)+"}^2}";
            }
            else{
                texOut<<"+c_{"+to_string(fillingIndex)+"}\\frac{N_"+to_string(speciesLabel)+"^2}{N_"+to_string(speciesLabel)+"^2+K_{"+to_string(fillingIndex)+"}^2}";
            }
            fillingIndex++;
        }
        texOut<<"-\\lambda_"+to_string(i)+"\\right)\\langle N_"+to_string(i)+"\\rangle";
        texOut<<"\\\\"<<endl;
    }
    texOut<<"\\end{aligned}"<<endl;
    texOut<<"\\end{equation}"<<endl;
    texOut.close();
    for(int i=0;i<numSpecies;i++){
        inData>>trueDecay[i];
    }
    inData.close();

    

    ofstream outPlot(inputFolder+"//trueOutPlot.txt");
    for(int i=0;i<(int)trueParameters.size();i++){
        for(int j=0;j<(int)trueParameters[i].size();j++){
            outPlot<<trueParameters[i][j]<<" ";
        }
    }
    for(int i=0;i<numSpecies;i++){
        outPlot<<trueDecay[i]<<" ";
    }
    outPlot.close();

    outPlot.open(inputFolder+"//fitData.txt");

    for(int fileIndex=0;fileIndex<numParticles;fileIndex++){
        inData.open(inputFolder+"//particleDump_globalBest_"+to_string(fileIndex)+".txt");
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
                inData>>trueParameters[fillingIndex][2];//power 
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