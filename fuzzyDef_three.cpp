#include <map>
#include <iostream>
#include <tuple>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>


#include "fuzzyDef_three.h"

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

Particle::Particle(vector<vector<hillStruct> > inSolution, tuple<double,double,double,double> bounds){
	sampleSolution=inSolution;
	int numSpecies=sampleSolution.size();
	normCurrentPos.resize(numSpecies);
	normBestPos.resize(numSpecies);
	normVelocity.resize(numSpecies);
	constCurrentPos.resize(numSpecies);
	constBestPos.resize(numSpecies);
	constVelocity.resize(numSpecies);
	powerCurrentPos.resize(numSpecies);
	powerBestPos.resize(numSpecies);
	powerVelocity.resize(numSpecies);
	decayConsts.resize(numSpecies);
	decayVelocities.resize(numSpecies);
	bestDecayConsts.resize(numSpecies);
	bestFitness=0;
	currentFitness=0;
	constBound=get<0>(bounds);
	normBound=get<1>(bounds);
    powerBound=get<2>(bounds);
	decayBound=get<3>(bounds);
	for(int i=0;i<(int)inSolution.size();i++){
		int numInteractions=inSolution[i].size();
		normCurrentPos[i].resize(numInteractions);
		normBestPos[i].resize(numInteractions);
		normVelocity[i].resize(numInteractions);
		constCurrentPos[i].resize(numInteractions);
		constBestPos[i].resize(numInteractions);
		constVelocity[i].resize(numInteractions);
		powerCurrentPos[i].resize(numInteractions);
		powerBestPos[i].resize(numInteractions);
		powerVelocity[i].resize(numInteractions);
	}
}

Particle::~Particle(){
	
	
	
}

void Particle::dumpParticleDetails(string rootFolder, string id){
	ofstream outParticle(rootFolder+"particleDump_"+id+".txt");
	outParticle<<"#m N: species Power Constant Norm Pos/Neg"<<endl;
	outParticle<<normBestPos.size()<<" ";
	int reactionCount(0);
	for(int i=0;i<(int)normBestPos.size();i++){
		for(int j=0;j<(int)normBestPos[i].size();j++){
			reactionCount++;
		}
	}
	outParticle<<reactionCount<<endl;
	
	for(int i=0;i<(int)sampleSolution.size();i++){
		outParticle<<i<<" "<<sampleSolution[i].size()<<" ";
		for(int j=0;j<(int)sampleSolution[i].size();j++){
			outParticle<<sampleSolution[i][j].speciesLabel<<" "<<sampleSolution[i][j].power<<" "<<sampleSolution[i][j].constant<<" "<<sampleSolution[i][j].normalization<<" "<<sampleSolution[i][j].hillBool<<" ";
		}
		outParticle<<endl;
	}
	
	for(int i=0;i<(int)decayConsts.size();i++){
		outParticle<<decayConsts[i]<<endl;
	}
	outParticle.close();
}



double Particle::performUpdate(boost::mt19937* inRand, FuzzyTree* fuzzyStruct, double parameterVectorToSend[]){
	double interDelta(0);
	
	int fillingIndex(0);
	for(int i=0;i<(int)normCurrentPos.size();i++){
		for(int j=0;j<(int)normCurrentPos[i].size();j++){
			double rand1((double)(*inRand)()/(double)(*inRand).max());
			double rand2((double)(*inRand)()/(double)(*inRand).max());
			double proposedVelocity(0);
			proposedVelocity+=(*fuzzyStruct).inertia*normVelocity[i][j];
			proposedVelocity+=(*fuzzyStruct).social*rand1*(parameterVectorToSend[fillingIndex]-normCurrentPos[i][j]);
			proposedVelocity+=(*fuzzyStruct).cognitive*rand2*(normBestPos[i][j]-normCurrentPos[i][j]);
			interDelta+=pow((parameterVectorToSend[fillingIndex]-normCurrentPos[i][j])/normBound,2);
			//Checks velocity limits
			if(proposedVelocity<(*fuzzyStruct).U*(-1.*normBound)){
				proposedVelocity=(*fuzzyStruct).U*(-1.*normBound);
			}
			if(proposedVelocity>(*fuzzyStruct).U*(normBound)){
				proposedVelocity=(*fuzzyStruct).U*(normBound);
			}
			//Checks bounds of parameter space;
			if(normCurrentPos[i][j]+proposedVelocity<0){
				normCurrentPos[i][j]=(-1.*proposedVelocity)*.1*((double)(*inRand)()/(double)(*inRand).max());
				normVelocity[i][j]=normCurrentPos[i][j];
			}
			else{
				if(normCurrentPos[i][j]+proposedVelocity>normBound){
					normCurrentPos[i][j]=normBound;
					normVelocity[i][j]=(-1.*proposedVelocity)*.1;
				}
				else{
					normCurrentPos[i][j]+=proposedVelocity;
					normVelocity[i][j]=proposedVelocity;
				}
			}
			fillingIndex++;
		}
	}
	
	
	for(int i=0;i<(int)constCurrentPos.size();i++){
		for(int j=0;j<(int)constCurrentPos[i].size();j++){
			double rand1((double)(*inRand)()/(double)(*inRand).max());
			double rand2((double)(*inRand)()/(double)(*inRand).max());
			double proposedVelocity(0);
			proposedVelocity+=(*fuzzyStruct).inertia*constVelocity[i][j];
			proposedVelocity+=(*fuzzyStruct).social*rand1*(parameterVectorToSend[fillingIndex]-constCurrentPos[i][j]);
			proposedVelocity+=(*fuzzyStruct).cognitive*rand2*(constBestPos[i][j]-constCurrentPos[i][j]);
			interDelta+=pow((parameterVectorToSend[fillingIndex]-constCurrentPos[i][j])/constBound,2);
			//Checks velocity limits
			if(proposedVelocity<(*fuzzyStruct).U*(-1.*constBound)){
				proposedVelocity=-1.*(*fuzzyStruct).U*(constBound);
			}
			if(proposedVelocity>(*fuzzyStruct).U*(constBound)){
				proposedVelocity=(*fuzzyStruct).U*(constBound);
			}
			//Checks bounds of parameter space;
			if(constCurrentPos[i][j]+proposedVelocity<0){
				constCurrentPos[i][j]=(-1.*proposedVelocity)*.1*((double)(*inRand)()/(double)(*inRand).max());
				constVelocity[i][j]=constCurrentPos[i][j];
			}
			else{
				if(constCurrentPos[i][j]+proposedVelocity>constBound){
					constCurrentPos[i][j]=constBound;
					constVelocity[i][j]=(-1.*proposedVelocity)*.1;
				}
				else{
					constCurrentPos[i][j]+=proposedVelocity;
					constVelocity[i][j]=proposedVelocity;
				}
			}
			fillingIndex++;
		}
	}

    for(int i=0;i<(int)powerCurrentPos.size();i++){
		for(int j=0;j<(int)powerCurrentPos[i].size();j++){
			double rand1((double)(*inRand)()/(double)(*inRand).max());
			double rand2((double)(*inRand)()/(double)(*inRand).max());
			double proposedVelocity(0);
			proposedVelocity+=(*fuzzyStruct).inertia*powerVelocity[i][j];
			proposedVelocity+=(*fuzzyStruct).social*rand1*(parameterVectorToSend[fillingIndex]-powerCurrentPos[i][j]);
			proposedVelocity+=(*fuzzyStruct).cognitive*rand2*(powerBestPos[i][j]-powerCurrentPos[i][j]);
			interDelta+=pow((parameterVectorToSend[fillingIndex]-powerCurrentPos[i][j])/powerBound,2);
			//Checks velocity limits
			if(proposedVelocity<(*fuzzyStruct).U*(-1.*powerBound)){
				proposedVelocity=-1.*(*fuzzyStruct).U*(powerBound);
			}
			if(proposedVelocity>(*fuzzyStruct).U*(powerBound)){
				proposedVelocity=(*fuzzyStruct).U*(powerBound);
			}
			//Checks bounds of parameter space;
			if(powerCurrentPos[i][j]+proposedVelocity<0){
				powerCurrentPos[i][j]=powerBound*.05;
				powerVelocity[i][j]=(-1.*proposedVelocity)*.1;
			}
			else{
				if(powerCurrentPos[i][j]+proposedVelocity>powerBound){
					powerCurrentPos[i][j]=powerBound;
					powerVelocity[i][j]=(-1.*proposedVelocity)*.1;
				}
				else{
					powerCurrentPos[i][j]+=proposedVelocity;
					powerVelocity[i][j]=proposedVelocity;
				}
			}
			fillingIndex++;
		}
	}
	
	for(int i=0;i<(int)decayConsts.size();i++){
		double rand1((double)(*inRand)()/(double)(*inRand).max());
		double rand2((double)(*inRand)()/(double)(*inRand).max());
		double proposedVelocity(0);
		proposedVelocity+=(*fuzzyStruct).inertia*decayVelocities[i];
		proposedVelocity+=(*fuzzyStruct).social*rand1*(parameterVectorToSend[fillingIndex]-decayConsts[i]);
		proposedVelocity+=(*fuzzyStruct).cognitive*rand2*(bestDecayConsts[i]-decayConsts[i]);
		interDelta+=pow((parameterVectorToSend[fillingIndex]-decayConsts[i])/decayBound,2);
		//Checks velocity limits
		if(proposedVelocity<(*fuzzyStruct).U*(-1.*decayBound)){
			proposedVelocity=-1.*(*fuzzyStruct).U*(decayBound);
		}
		if(proposedVelocity>(*fuzzyStruct).U*(decayBound)){
			proposedVelocity=(*fuzzyStruct).U*(decayBound);
		}
		//Checks bounds of parameter space;
		if(decayConsts[i]+proposedVelocity<0){
			decayConsts[i]+=(-1.*proposedVelocity);
			decayVelocities[i]=(-1.*proposedVelocity)*.25;
		}
		else{
			if(decayConsts[i]+proposedVelocity>normBound){
				decayConsts[i]=normBound-proposedVelocity;
				decayVelocities[i]=(-1.*proposedVelocity)*.25;
			}
			else{
				decayConsts[i]+=proposedVelocity;
				decayVelocities[i]=proposedVelocity;
			}
		}
		fillingIndex++;
	}
	
	return interDelta;
}