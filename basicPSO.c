/* Very basic particle swarm optimization setup
    Simulates the following equations:
    dA/dt=-L_3*[C]
    dB/dt=L_1*[A]*[B]
    dC/dt=L_2*[B]*[C]-L_4*C

    The parameters to optimize are L_1, L_2, and L_3. L_4 is inserted to give the system a steady state to reach.

*/
//Warnings, using rand() from the stdio library. This is terrible from the point of view of the stochastics, so any serious reimplmentation should use a MT19937 or better algorithm to produce the random numbers. 

//Compiles with gcc, there are no external files, and all includes are part of the std library.




#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct _PSO{

    double currentFitness, bestFitness;
    double curAPromo, curBPromo, curCInhib;
    double curAVel, curBVel, curCVel;
    double L, U, cCog, cSoc, w;
    double ABound, BBound, CBound;
    double cDecayConstant;
    double* particleBestSet;

} PSOStruct;

double sgn(double in){
    return (in<0) ? -1 : ((in>0)? 1 : 0);
}

void updatePSOPosition(PSOStruct* PSOObject, double* globalBestSet, double* particleBestSet);
void rungeKuttaIteration(double* speciesIn, double deltaT, double currentTime, double stoppingTime, const int speciesSize, PSOStruct* currentParticle);
double calculateFitnessFunction(double* testIn, double* fitIn, const int dataSize);
void testSolutions(double* speciesIn, double deltaT, PSOStruct* testParticle, PSOStruct* fitParticle);

int main(){
    srand(time(NULL));
    const int numOfParameters=3;
    double globalBestSet[3]={2,1,0.45};
    //Creates bounds for the parameters
    double ABound, BBound, CBound;
    ABound=1;
    BBound=1;
    CBound=1;
    double CDecayConstant=0.5;

    //Initializes array lengths and particle numbers
    const int numParticles=10;
    const int speciesNumber=3;
    PSOStruct particleArray[numParticles];

    //Initialize test data
    PSOStruct testParticle;
    testParticle.curAPromo=.2;
    testParticle.curBPromo=.1;
    testParticle.curCInhib=0.5;
    testParticle.cDecayConstant=CDecayConstant;

    //Creates the data to train the particle swarm and stores it in testArray
    double stoppingTimes[6]={0,2,4,8,16,32};
    const int numberOfReports=5;
    double timeIncrement=0.0002;
    double testArray[5][3];
    for(int i=0;i<5;i++){
        for(int j=0;j<3;j++){
            testArray[i][j]=0;
        }
    }
    double testSpecies[3]={2,1,.45};
    for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
        rungeKuttaIteration(testSpecies, timeIncrement, stoppingTimes[reportIndex], stoppingTimes[reportIndex+1], speciesNumber, &testParticle);
        for(int i=0;i<speciesNumber;i++){
            testArray[reportIndex][i]=testSpecies[i];
        }
    }
    


    for(int particleInit=0;particleInit<numParticles;particleInit++){
        //Random initialization for each particle, bad statistics on assignment. Should be improved in a more sophisticated implementation.
        double initCurFitness=1e24;
        double initBestFitness=0;
        double initCurAPromo=(double)rand()/RAND_MAX*ABound;
        double initCurBPromo=(double)rand()/RAND_MAX*BBound;
        double initCurCInhib=(double)rand()/RAND_MAX*CBound;
        double initCurAVel=0;
        double initCurBVel=0;
        double initCurCVel=0;
        double initL=0;
        double initU=.2;
        double initcCog=2;
        double initcSoc=2;
        double initW=.6;
        double initABound=ABound;
        double initBBound=BBound;
        double initCBound=CBound;
        double cDecayConstant=CDecayConstant;
        double initBestSet[3]={initCurAPromo,initCurBPromo,initCurCInhib};

        PSOStruct interParticle={initCurFitness,initBestFitness,initCurAPromo,initCurBPromo,initCurCInhib,initCurAVel,initCurBVel,initCurCVel,initL,initU,initcCog,initcSoc,initW,initABound,initBBound,initCBound,cDecayConstant,initBestSet};

        particleArray[particleInit]=interParticle;
    }

    double globalBestFitness=0;
    double bestFitnessValues[numParticles];
    for(int i=0;i<numParticles;i++){
        bestFitnessValues[i]=0;
    }

    
    //Initialize particles and calculate their initial fitness values and positions in the space
    for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
        PSOStruct* interParticle=&particleArray[currentParticle];
        double speciesArray[3]={2,1,.45};
        double errorValue=0;
        for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
            rungeKuttaIteration(speciesArray, timeIncrement, stoppingTimes[reportIndex], stoppingTimes[reportIndex+1], speciesNumber, interParticle);
            errorValue+=calculateFitnessFunction(testArray[reportIndex],speciesArray,3);
        }
        (*interParticle).bestFitness=errorValue;
        (*interParticle).currentFitness=errorValue;
        if(errorValue<(*interParticle).bestFitness){
            (*interParticle).bestFitness=errorValue;
            (*interParticle).particleBestSet[0]=(*interParticle).curAPromo;
            (*interParticle).particleBestSet[1]=(*interParticle).curBPromo;
            (*interParticle).particleBestSet[2]=(*interParticle).curCInhib;
        }
        bestFitnessValues[currentParticle]=(*interParticle).bestFitness;
    }
    
    //Find the global best position amongst the particles
    int bestParticleIndex=0;
    double compareHold=1e10;
    for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
        if(compareHold>bestFitnessValues[currentParticle]){
            compareHold=bestFitnessValues[currentParticle];
            bestParticleIndex=currentParticle;
        }
    }

    //Set the best parameter set in array to share amongst the particles
    globalBestSet[0]=particleArray[bestParticleIndex].curAPromo;
    globalBestSet[1]=particleArray[bestParticleIndex].curBPromo;
    globalBestSet[2]=particleArray[bestParticleIndex].curCInhib;

    //Performs PSO update for each particle
    for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
        updatePSOPosition(&particleArray[currentParticle], globalBestSet, particleArray[currentParticle].particleBestSet);
    }

    //Iterate PSO using the same steps as above
    const int numIterations=100;
    for(int iteration=0;iteration<numIterations;iteration++){
        for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
            PSOStruct* interParticle=&particleArray[currentParticle];
            double speciesArray[3]={2,1,.45};
            double errorValue=0;
            for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
                rungeKuttaIteration(speciesArray, timeIncrement, stoppingTimes[reportIndex], stoppingTimes[reportIndex+1], speciesNumber, interParticle);
                errorValue+=calculateFitnessFunction(testArray[reportIndex],speciesArray,3);
            }
            (*interParticle).currentFitness=errorValue;
            printf("%f ",errorValue);
            if(errorValue<(*interParticle).bestFitness){
                (*interParticle).bestFitness=errorValue;
                (*interParticle).particleBestSet[0]=(*interParticle).curAPromo;
                (*interParticle).particleBestSet[1]=(*interParticle).curBPromo;
                (*interParticle).particleBestSet[2]=(*interParticle).curCInhib;
            }
            bestFitnessValues[currentParticle]=(*interParticle).bestFitness;
        }
            
        int bestParticleIndex=0;
        double compareHold=1e10;
        for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
            if(compareHold>bestFitnessValues[currentParticle]){
                compareHold=bestFitnessValues[currentParticle];
                bestParticleIndex=currentParticle;
            }
        }

        globalBestSet[0]=particleArray[bestParticleIndex].curAPromo;
        globalBestSet[1]=particleArray[bestParticleIndex].curBPromo;
        globalBestSet[2]=particleArray[bestParticleIndex].curCInhib;

        for(int currentParticle=0;currentParticle<numParticles;currentParticle++){
            updatePSOPosition(&particleArray[currentParticle], globalBestSet, particleArray[currentParticle].particleBestSet);
        }
        printf("%f\n",compareHold);
    }

    //Print out results to both console and to files
    PSOStruct bestParticle=particleArray[0];
    bestParticle.curAPromo=globalBestSet[0];
    bestParticle.curBPromo=globalBestSet[1];
    bestParticle.curCInhib=globalBestSet[2];
    printf("%f %f %f\n",testParticle.curAPromo,testParticle.curBPromo,testParticle.curCInhib);
    printf("%f %f %f\n",globalBestSet[0],globalBestSet[1],globalBestSet[2]);
    double solutionsTest[3]={2,1,.45};
    testSolutions(solutionsTest,timeIncrement,&testParticle,&bestParticle);

    return 0;
}

void updatePSOPosition(PSOStruct* PSOObject, double* globalBestSet, double* particleBestSet){

    double proposedUpdate=0;
    //Performs A update
    proposedUpdate+=(*PSOObject).w*(*PSOObject).curAVel;
    proposedUpdate+=(*PSOObject).cSoc*rand()/((double)RAND_MAX)*(globalBestSet[0]-(*PSOObject).curAPromo);
    proposedUpdate+=(*PSOObject).cCog*rand()/((double)RAND_MAX)*(particleBestSet[0]-(*PSOObject).curAPromo);
    
    //Checks velocity magnitude
    if(fabs(proposedUpdate)>(*PSOObject).U*(*PSOObject).ABound){
        proposedUpdate=sgn(proposedUpdate)*(*PSOObject).U*(*PSOObject).ABound;
    }

    //Checks positive bounds
    if((proposedUpdate+(*PSOObject).curAPromo)>(*PSOObject).ABound){
        (*PSOObject).curAPromo=(*PSOObject).ABound;
        (*PSOObject).curAVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
    } else{
        if((proposedUpdate+(*PSOObject).curAPromo)<0){
            (*PSOObject).curAPromo=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
            (*PSOObject).curAVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
        } else{
            (*PSOObject).curAPromo+=proposedUpdate;
            (*PSOObject).curAVel=proposedUpdate;
        }
    }

    proposedUpdate=0;

    //Performs B update
    proposedUpdate+=(*PSOObject).w*(*PSOObject).curBVel;
    proposedUpdate+=(*PSOObject).cSoc*rand()/((double)RAND_MAX)*(globalBestSet[1]-(*PSOObject).curBPromo);
    proposedUpdate+=(*PSOObject).cCog*rand()/((double)RAND_MAX)*(particleBestSet[1]-(*PSOObject).curBPromo);
    
    //Checks velocity magnitude
    if(fabs(proposedUpdate)>(*PSOObject).U*(*PSOObject).BBound){
        proposedUpdate=sgn(proposedUpdate)*(*PSOObject).U*(*PSOObject).BBound;
    }

    //Checks positive bounds
    if((proposedUpdate+(*PSOObject).curBPromo)>(*PSOObject).BBound){
        (*PSOObject).curBPromo=(*PSOObject).BBound;
        (*PSOObject).curBVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
    } else{
        if((proposedUpdate+(*PSOObject).curBPromo)<0){
            (*PSOObject).curBPromo=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
            (*PSOObject).curBVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
        } else{
            (*PSOObject).curBPromo+=proposedUpdate;
            (*PSOObject).curBVel=proposedUpdate;
        }
    }

    proposedUpdate=0;

    //Performs C update
    proposedUpdate+=(*PSOObject).w*(*PSOObject).curCVel;
    proposedUpdate+=(*PSOObject).cSoc*rand()/((double)RAND_MAX)*(globalBestSet[2]-(*PSOObject).curCInhib);
    proposedUpdate+=(*PSOObject).cCog*rand()/((double)RAND_MAX)*(particleBestSet[2]-(*PSOObject).curCInhib);
    
    //Checks velocity magnitude
    if(fabs(proposedUpdate)>(*PSOObject).U*(*PSOObject).CBound){
        proposedUpdate=sgn(proposedUpdate)*(*PSOObject).U*(*PSOObject).CBound;
    }

    //Checks positive bounds
    if((proposedUpdate+(*PSOObject).curCInhib)>(*PSOObject).CBound){
        (*PSOObject).curCInhib=(*PSOObject).BBound;
        (*PSOObject).curCVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
    } else{
        if((proposedUpdate+(*PSOObject).curCInhib)<0){
            (*PSOObject).curCInhib=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
            (*PSOObject).curCVel=-1.*sgn(proposedUpdate)*proposedUpdate/10.;
        } else{
            (*PSOObject).curCInhib+=proposedUpdate;
            (*PSOObject).curCVel=proposedUpdate;
        }
    }
}

//Calculates mean squared difference in the data at given time point
double calculateFitnessFunction(double* testIn, double* fitIn, const int dataSize){
    double outHold=0;
    for(int i=0;i<dataSize;i++){
        outHold+=pow(testIn[i]-fitIn[i],2);
    }
    return outHold;
}

//Performs fourth order runge kutta iteration from currentTime to stoppingTime in increments of deltaT using parameters taken from currentParticle and places the results in array speciesIn
void rungeKuttaIteration(double* speciesIn, double deltaT, double currentTime, double stoppingTime, const int speciesSize, PSOStruct* currentParticle){
	int n=(int)(((stoppingTime-currentTime))/(deltaT));
    static double species[3];
	for(int t=0;t<n;t++){
		double k1[speciesSize];
		k1[0]=deltaT*(-1.*(*currentParticle).curCInhib)*speciesIn[2];
        k1[1]=deltaT*((*currentParticle).curAPromo)*speciesIn[0]*speciesIn[1];
        k1[2]=deltaT*((*currentParticle).curBPromo*speciesIn[1]*speciesIn[2]-(*currentParticle).cDecayConstant*speciesIn[2]);
		double k2[speciesSize];
		k2[0]=deltaT*(-1.*(*currentParticle).curCInhib)*(speciesIn[2]+k1[2]/2.);
        k2[1]=deltaT*((*currentParticle).curAPromo)*(speciesIn[0]+k1[0]/2.)*(speciesIn[1]+k1[1]/2.);
        k2[2]=deltaT*((*currentParticle).curBPromo*(speciesIn[1]+k1[1]/2.)*(speciesIn[2]+k1[2]/2.)-(*currentParticle).cDecayConstant*(speciesIn[2]+k1[2]/2.));
		double k3[speciesSize];
		k3[0]=deltaT*(-1.*(*currentParticle).curCInhib)*(speciesIn[2]+k2[2]/2.);
        k3[1]=deltaT*((*currentParticle).curAPromo)*(speciesIn[0]+k2[0]/2.)*(speciesIn[1]+k2[1]/2.);
        k3[2]=deltaT*((*currentParticle).curBPromo*(speciesIn[1]+k2[1]/2.)*(speciesIn[2]+k2[2]/2.)-(*currentParticle).cDecayConstant*(speciesIn[2]+k2[2]/2.));
		double k4[speciesSize];
		k4[0]=deltaT*(-1.*(*currentParticle).curCInhib)*(speciesIn[2]+k3[2]);
        k4[1]=deltaT*((*currentParticle).curAPromo)*(speciesIn[0]+k3[0])*(speciesIn[1]+k3[1]);
        k4[2]=deltaT*((*currentParticle).curBPromo)*(speciesIn[1]+k3[1])*(speciesIn[2]+k3[2])-deltaT*(*currentParticle).cDecayConstant*(speciesIn[2]+k3[2]);
		for(int i=0;i<speciesSize;i++){
			speciesIn[i]+=k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.;
		}
	}
}

//Outputs the testing data to "testOut.txt" by iterating over the test particle parameters using rungeKutta, and output the fit data into "fitOut.txt" using the fit particle parameters.
void testSolutions(double* speciesIn, double deltaT, PSOStruct* testParticle, PSOStruct* fitParticle){
    FILE *fp;
    fp=fopen("testOut.txt","w");
    double stoppingTimes[6]={0,2,4,8,16,32};
    double speciesHold[3];
    for(int i=0;i<3;i++){
        speciesHold[i]=speciesIn[i];
    }
    const int numberOfReports=5;
    for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
        rungeKuttaIteration(speciesIn, deltaT, stoppingTimes[reportIndex], stoppingTimes[reportIndex+1], 3, testParticle);
        fprintf(fp,"%f ",stoppingTimes[reportIndex+1]);
        for(int i=0;i<3;i++){
            fprintf(fp,"%f ",speciesIn[i]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    fp=fopen("fitOut.txt","w");
    for(int reportIndex=0;reportIndex<numberOfReports;reportIndex++){
        rungeKuttaIteration(speciesHold, deltaT, stoppingTimes[reportIndex], stoppingTimes[reportIndex+1], 3, fitParticle);
        fprintf(fp,"%f ",stoppingTimes[reportIndex+1]);
        for(int i=0;i<3;i++){
            fprintf(fp,"%f ",speciesHold[i]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

