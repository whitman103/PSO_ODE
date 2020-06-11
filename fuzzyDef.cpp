#include <map>
#include <tuple>
#include <string>
#include <math.h>
#include <vector>


#include "fuzzyDef.h"

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
	
}

	