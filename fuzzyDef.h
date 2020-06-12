#ifndef FUZZY_H
#define FUZZY_H

#include <map>
#include <tuple>
#include <string>
#include <vector>

using namespace std;

class FuzzyTree{
	public:
	FuzzyTree(double inDelta);
	~FuzzyTree();
	
	double phi;
	double delta, delta1, delta2, delta3, deltaMax;
	//phi Membership values are set as better, same, worse
	vector<double> phiMembershipValues;
	
	
	//delta Membership values are set as same, near, far
	vector<double> deltaMembershipValues;
	
	double inertia;
	double social;
	double cognitive;
	double L;
	double U;
	map<string,double> inertiaMap;
	map<string,double> socialMap;
	map<string,double> cognitiveMap;
	map<string,double> LMap;
	map<string,double> UMap;
	map<int,string> linguisticMap;
	
	void setParameters();
	void setInertia();
	void setSocial();
	void setCognitive();
	void setL();
	void setU();
	void calculatePhiMembershipValues();
	void calculateDeltaMembershipValues();
	
	
	
	private:
	
};

typedef struct{
	bool hillBool;
	double constant;
	double power;
	double normalization;
	int speciesLabel;
} hillStruct;

class Particle{
	public:
	Particle(vector<vector<hillStruct> > inSolution, tuple<double,double,double> bounds);
	~Particle();
	vector<vector<hillStruct> > sampleSolution;
	vector<vector<double> > normCurrentPos;
	vector<vector<double> > normBestPos;
	vector<vector<double> > normVelocity;
	vector<vector<double> > constCurrentPos;
	vector<vector<double> > constBestPos;
	vector<vector<double> > constVelocity;
	vector<double> decayConsts;
	vector<double> decayVelocities;
	vector<double> bestDecayConsts;
	double currentFitness;
	double lastFitness;
	tuple<double,double> constBounds;
	tuple<double,double> normBounds;
	double decayBound;
	private:
};



#endif