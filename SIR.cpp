#include <map>
#include <iostream>
#include <tuple>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
using namespace std;


#include "SIR.h"

//Values taken directly from Nobile2017 paper
FuzzyTree::FuzzyTree(double inDelta){
	inertiaMap["Low"]=0.3;
	inertiaMap["Medium"]=0.5;
	inertiaMap["High"]=1.0;
	socialMap["Low"]=1.0;
	socialMap["Medium"]=2.0;
	socialMap["High"]=3.0;
	cognitiveMap["Low"]=0.1;
	cognitiveMap["Medium"]=1.5;
	cognitiveMap["High"]=3.0;
	LMap["Low"]=0.0;
	LMap["Medium"]=0.001;
	LMap["High"]=0.01;
	UMap["Low"]=0.1;
	UMap["Medium"]=0.15;
	UMap["High"]=0.2;
	linguisticMap[0]="Low";
	linguisticMap[1]="Medium";
	linguisticMap[2]="High";
	
	deltaMax=inDelta;
	delta1=0.2*deltaMax;
	delta2=0.4*deltaMax;
	delta3=0.6*deltaMax;
	
	inertia=inertiaMap["Medium"];
	social=socialMap["Medium"];
	cognitive=cognitiveMap["Medium"];
	L=LMap["Medium"];
	U=UMap["Medium"];
};


FuzzyTree::~FuzzyTree(){
};

void FuzzyTree::setParameters(){
	
	calculateDeltaMembershipValues();
	calculatePhiMembershipValues();
	setInertia();
	setSocial();
	setCognitive();
	setL();
	setU();
}

void FuzzyTree::calculateDeltaMembershipValues(){
	if(delta<delta1){
		deltaMembershipValues={1,0,0};
	} 
	else{
		if(delta>delta1&&delta<delta2){
			deltaMembershipValues={(delta2-delta)/(delta2-delta1),(delta-delta1)/(delta2-delta1),0};
		}
		else{
			if(delta>delta2&&delta<delta3){
				deltaMembershipValues={0,(delta3-delta)/(delta3-delta2),(delta-delta2)/(delta3-delta2)};
			}
			else{
				if(delta>delta3){
					deltaMembershipValues={0,0,1};
				}
			}
		}
	}
}

void FuzzyTree::calculatePhiMembershipValues(){
	if(phi<0){
		phiMembershipValues={-phi,1-fabs(phi),0};
	}
	else{
		phiMembershipValues={0,1-fabs(phi),phi};
	}
}

void FuzzyTree::setInertia(){
	inertia=0;
	inertia+=inertiaMap["Low"]*(phiMembershipValues[2]+deltaMembershipValues[0]);
	inertia+=inertiaMap["Medium"]*(phiMembershipValues[1]+deltaMembershipValues[1]);
	inertia+=inertiaMap["High"]*(phiMembershipValues[0]+deltaMembershipValues[2]);
	inertia/=2.;
}

void FuzzyTree::setSocial(){
	social=0;
	social+=socialMap["Low"]*(phiMembershipValues[0]+deltaMembershipValues[1]);
	social+=socialMap["Medium"]*(phiMembershipValues[1]+deltaMembershipValues[0]);
	social+=socialMap["High"]*(phiMembershipValues[2]+deltaMembershipValues[2]);
	social/=2.;
}

void FuzzyTree::setCognitive(){
	cognitive=0;
	cognitive+=cognitiveMap["Low"]*deltaMembershipValues[2];
	cognitive+=cognitiveMap["Medium"]*(phiMembershipValues[2]+phiMembershipValues[1]+deltaMembershipValues[0]+deltaMembershipValues[1]);
	cognitive+=cognitiveMap["High"]*(phiMembershipValues[0]);
	cognitive/=2.;
}

void FuzzyTree::setL(){
	L=0;
	L+=LMap["Low"]*(phiMembershipValues[1]+phiMembershipValues[0]+deltaMembershipValues[2]);
	L+=LMap["Medium"]*(deltaMembershipValues[0]+deltaMembershipValues[1]);
	L+=LMap["High"]*(phiMembershipValues[2]);
	L/=2.;
}

void FuzzyTree::setU(){
	U=0;
	U+=UMap["Low"]*(deltaMembershipValues[0]);
	U+=UMap["Medium"]*(phiMembershipValues[0]+phiMembershipValues[1]+deltaMembershipValues[1]);
	U+=UMap["High"]*(phiMembershipValues[2]+deltaMembershipValues[2]);
	U/=2.;
}

void FuzzyTree::calculatePhi(double lastFitness, double currentFitness){
	phi=0;
	phi=delta/deltaMax;
	phi*=(min(currentFitness,phiNormalization)-min(lastFitness,phiNormalization))/phiNormalization;
}

Particle::Particle(int numOfParameters, vector<double> initBounds, vector<double (*)(Particle*,vector<double>&)> initFunctions){
    currentSolution.resize(numOfParameters);
    bestSolution.resize(numOfParameters);
    currentVelocity.resize(numOfParameters);
    interactionFunctions=initFunctions;
    bestFitness=0;
    currentFitness=0;
	
}

Particle::~Particle(){
}

void Particle::dumpParticleDetails(string rootFolder, string id){
}

double sgn(double in){
    return (in<0) ? -1 : ((in>0)? 1 : 0);
}



double Particle::performUpdate(boost::mt19937* inRand, vector<double>& globalBest, FuzzyTree* fuzzyStruct){
    for(int i=0;i<(int)currentSolution.size();i++){
        double proposedUpdate(0);
        double rand1((double)(*inRand)()/(double)(*inRand).max());
        double rand2((double)(*inRand)()/(double)(*inRand).max());
        proposedUpdate+=(*fuzzyStruct).inertia*currentVelocity[i];
        proposedUpdate+=(*fuzzyStruct).social*rand1*(globalBest[i]-currentSolution[i]);
        proposedUpdate+=(*fuzzyStruct).cognitive*rand2*(bestSolution[i]-currentSolution[i]);
        if(fabs(proposedUpdate)>(*fuzzyStruct).U*bounds[i]){
            proposedUpdate=(*fuzzyStruct).U*bounds[i]*sgn(proposedUpdate);
        }

        if((currentSolution[i]+proposedUpdate)<0){
            currentSolution[i]=1./100.*bounds[i];
            currentVelocity[i]=-1.*proposedUpdate*1./10.;
        } else{
            if((currentSolution[i]+proposedUpdate)>bounds[i]){
                currentSolution[i]=bounds[i];
                currentVelocity[i]=-1.*proposedUpdate*1./10.;
            } else{
                currentSolution[i]+=proposedUpdate;
                currentVelocity[i]=proposedUpdate;
            }
        }
    }
	return 0;
}

void Particle::unwrapParameters(){
    beta=currentSolution[0];
    delta=currentSolution[1];
    c=currentSolution[2];
    p=currentSolution[3];
}

//Species are currently T, I, V

double firstInteraction(Particle* currentParticle, vector<double>& species){
    return -1.*(*currentParticle).beta*species[0]*species[1];
}

double secondInteraction(Particle* currentParticle, vector<double>& species){
    return (*currentParticle).beta*species[0]*species[1]-(*currentParticle).c*species[2];
}

double thirdInteraction(Particle* currentParticle, vector<double>& species){
    return (*currentParticle).p*species[1]-(*currentParticle).delta*species[2];
}

void rungeKuttaUpdate(Particle* currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT){
    int numSpecies=speciesVec.size();
    (*currentParticle).unwrapParameters();
    int n=(int)(((stoppingTime-currentTime))/(deltaT));

    vector<double (*)(Particle*,vector<double>&)> interactionPointer=(*currentParticle).interactionFunctions;

    for(int t=0;t<n;t++){
        vector<double> interSpecies(numSpecies,0);
        vector<double> k1(numSpecies,0);
        for(int i=0;i<(int)k1.size();i++){
            k1[i]=interactionPointer[i](currentParticle,speciesVec);
        }
        transform(speciesVec.begin(),speciesVec.end(),interSpecies.begin(),interSpecies.begin(),std::plus<double>());
        vector<double> k2(numSpecies,0);
        for(int i=0;i<(int)k2.size();i++){
            k2[i]=interactionPointer[i](currentParticle,interSpecies)/2.;
        }
        interSpecies=speciesVec;
        transform(speciesVec.begin(),speciesVec.end(),interSpecies.begin(),interSpecies.begin(),std::plus<double>());
        vector<double> k3(numSpecies,0);
        for(int i=0;i<(int)k3.size();i++){
            k3[i]=interactionPointer[i](currentParticle,interSpecies)/2.;
        }
        interSpecies=speciesVec;
        transform(speciesVec.begin(),speciesVec.end(),interSpecies.begin(),interSpecies.begin(),std::plus<double>());
        vector<double> k4(numSpecies,0);
        for(int i=0;i<(int)k4.size();i++){
            k4[i]=interactionPointer[i](currentParticle,interSpecies);
        }
        for(int i=0;i<numSpecies;i++){
            speciesVec[i]+=k1[i]/6.+k2[i]*2./3.+k3[i]*2./3.+k4[i]*6.;
        }
    }
}

double fitnessFunction(vector<vector<double> >& trueIn, vector<vector<double> >& testIn){
    double interHold(0);
    for(int i=0;i<(int)trueIn.size();i++){
        for(int j=0;j<(int)trueIn[i].size();j++){
            interHold+=pow(trueIn[i][j]-testIn[i][j],2);
        }
    }
    return interHold;

}