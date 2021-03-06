#ifndef SIR_H
#define SIR_H

#include <map>
#include <tuple>
#include <string>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
using boost::mt19937;



class Particle{
    public:
    Particle(int numOfParameters, vector<double> initBounds, vector<double (*)(Particle*,vector<double>&)> initFunctions);
    vector<double> currentSolution;
    vector<double> bestSolution;
    vector<double> currentVelocity;
    vector<double> bounds;
    vector<double (*)(Particle*,vector<double>&)> interactionFunctions;
    double bestFitness;
    double currentFitness;
    double beta, delta, c, p;
    void unwrapParameters();

    void dumpParticleDetails(string rootFolder, string id);
    double performUpdate(boost::mt19937* inRand, vector<double>& globalBest, FuzzyTree* fuzzyStruct);

    private:
};

double sgn(double in);

double firstInteraction(Particle* currentParticle, vector<double>& species);

double secondInteraction(Particle* currentParticle, vector<double>& species);

double thirdInteraction(Particle* currentParticle, vector<double>& species);

void rungeKuttaUpdate(Particle* currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT);

double fitnessFunction(vector<vector<double> >& trueIn, vector<vector<double> >& testIn);

class FuzzyTree{
	public:
	FuzzyTree(double inDelta);
	~FuzzyTree();
	
	double phi, phiNormalization;
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
	
	void calculatePhi(double lastFitness, double currentFitness);
	
	
	
	private:
	
};

#endif