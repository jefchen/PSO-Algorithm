/*
 * main.cpp
 *
 *  Created on: Jul 19, 2012
 *      Author: jeffreychen_55
 */

using namespace std;
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <float.h>

const float WEIGHT = 0.792;
const float C_PERSONAL = 2;
const float C_NEIGHBOR = 2;
const bool USE_CONSTRICTION = true;

const int SWARM_SIZE = 30;
const int NEIGHBORHOOD_SIZE = 3;

const int ITERATION_MIN = 10;
const int ITERATION_MAX = 10000;
const int RUN_COUNT = 10;

// problem 1
/*
const int VAR_COUNT = 2;
const int RANGE[VAR_COUNT] = {3, 2};
float evaluate(float *x) {
	float x2 = pow(x[0], 2);
	float x4 = pow(x[0], 4);
	float y2 = pow(x[1], 2);
	float result = (4 - 2.1*x2 + x4/3)*x2 + x[0]*x[1] + (-4 + 4*y2)*y2;
	return result;
}
const float VMAX[VAR_COUNT] = {FLT_MAX, FLT_MAX};
*/

// problem 2
const int VAR_COUNT = 10;
const int RANGE[VAR_COUNT] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
float evaluate(float *problem) {
	float result = 0;
	for (int i=0; i<VAR_COUNT; i++) {
		result += pow(problem[i], 2);
	}
	return abs(result);
}
const float VMAX[VAR_COUNT] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
//const float VMAX[VAR_COUNT] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};


typedef struct Particle {
	float *x;
	float *v;
	float *pBest;
	float pBestResult;
	float *nBest;
	float nBestResult;
} Particle;

// copy
void copyArray(float *x, float *newX) {
	for (int j=0; j<VAR_COUNT; j++) {
		newX[j] = x[j];
	}
}

// generate random number between 0 and 1
float randomZeroAndOne()
{
	float scale=RAND_MAX+1.;
	float base=rand()/scale;
	float fine=rand()/scale;
	return base+fine/scale;
}

// problem 1
void initialize(Particle particles[SWARM_SIZE]) {
	float nBest[VAR_COUNT];
	float nBestResult = FLT_MAX;

	for (int i=0; i<SWARM_SIZE; i++) {
		// initialize memory
		particles[i].x = new float[VAR_COUNT];
		particles[i].v = new float[VAR_COUNT];
		particles[i].nBest = new float[VAR_COUNT];
		particles[i].pBest = new float[VAR_COUNT];

		for (int j=0; j<VAR_COUNT; j++) {
			// initialize position randomly
			particles[i].x[j] = (float)(rand() % (2000*RANGE[j]))/(float)1000 - 1.0*RANGE[j];

			// initialize velocity to 0
			particles[i].v[j] = 0;
		}

		// set current Best
		particles[i].pBestResult = evaluate(particles[i].x);
		copyArray(particles[i].x, particles[i].pBest);
	}

	// update nBest for all particles
	for (int i=0; i<SWARM_SIZE; i++) {
		nBestResult = FLT_MAX;
		for (int j=i-NEIGHBORHOOD_SIZE; j<i+NEIGHBORHOOD_SIZE+1; j++) {
			if (particles[(j+SWARM_SIZE) % SWARM_SIZE].pBestResult < nBestResult) {
				nBestResult = particles[(i+j+SWARM_SIZE) % SWARM_SIZE].pBestResult;
				copyArray(particles[(i+j+SWARM_SIZE) % SWARM_SIZE].pBest, nBest);
			}
		}
		particles[i].nBestResult = nBestResult;
		copyArray(nBest, particles[i].nBest);
	}
}

void updateVelocity(int i, Particle particles[SWARM_SIZE]) {
	float r1;
	float r2;
	for (int j=0; j<VAR_COUNT; j++) {
		// set r1 and r2
		r1 = randomZeroAndOne();
		r2 = randomZeroAndOne();
		if (USE_CONSTRICTION) {
			float theta = C_PERSONAL + C_NEIGHBOR;
			float factor = abs(2-theta-sqrt(pow(theta,2)-4*theta));
			factor = 2/factor;
			particles[i].v[j] += C_PERSONAL*r1*(particles[i].pBest[j]-particles[i].x[j]);
			particles[i].v[j] += C_NEIGHBOR*r2*(particles[i].nBest[j]-particles[i].x[j]);
			particles[i].v[j] *= factor;
		} else {
			particles[i].v[j] = WEIGHT*particles[i].v[j];
			particles[i].v[j] += C_PERSONAL*r1*(particles[i].pBest[j]-particles[i].x[j]);
			particles[i].v[j] += C_NEIGHBOR*r2*(particles[i].nBest[j]-particles[i].x[j]);
		}
		if (abs(particles[i].v[j]) > VMAX[j]) {
			particles[i].v[j] = particles[i].v[j] > 0 ? VMAX[j] : -VMAX[j];
		}
	}
}

void updatePosition(int i, Particle particles[SWARM_SIZE]) {
	for (int j=0; j<VAR_COUNT; j++) {
		particles[i].x[j] += particles[i].v[j];
	}
}

int main(int argc, char* argv[]) {

	// initialize random seed
	srand(time(NULL));

	for (int run=0; run < RUN_COUNT; run++) {
		// print run
		cout << "***** RUN " << run << " *********" << endl;

		// initialize swarm for question 1
		Particle particles[SWARM_SIZE];
		initialize(particles);

		// while termination criteria is not met
		float nBestResult = FLT_MAX;
		float nBest[VAR_COUNT];
		float pBestTemp = FLT_MAX;
		float gBestResult = FLT_MAX;
		float gBest[VAR_COUNT];
		float average =0;
		int iteration = 0;
		int printCount = ITERATION_MIN;
		while (iteration < ITERATION_MAX) {
			// for each particle
			for (int i=0; i<SWARM_SIZE; i++) {
				// update velocity
				updateVelocity(i, particles);

				// update position
				updatePosition(i, particles);

				// update particle's personal best
				pBestTemp = evaluate(particles[i].x);
				if (pBestTemp < particles[i].pBestResult) {
					particles[i].pBestResult = pBestTemp;
					copyArray(particles[i].x, particles[i].pBest);
				}

				// update global best for printing
				if (particles[i].pBestResult < gBestResult) {
					gBestResult = particles[i].pBestResult;
					copyArray(particles[i].pBest, gBest);
				}
			}

			// update nBest for all particles
			for (int i=0; i<SWARM_SIZE; i++) {
				nBestResult = FLT_MAX;
				for (int j=i-NEIGHBORHOOD_SIZE; j<i+NEIGHBORHOOD_SIZE+1; j++) {
					if (particles[(j+SWARM_SIZE) % SWARM_SIZE].pBestResult < nBestResult) {
						nBestResult = particles[(i+j+SWARM_SIZE) % SWARM_SIZE].pBestResult;
						copyArray(particles[(i+j+SWARM_SIZE) % SWARM_SIZE].pBest, nBest);
					}
				}
				particles[i].nBestResult = nBestResult;
				copyArray(nBest, particles[i].nBest);
			}

			// increment count
			iteration++;

			// print result
			if (iteration == printCount) {
				cout << "ITERATION: " << iteration << ",ACHIEVED: " << gBestResult << ", SOLUTION: [";
				for (int i=0; i<VAR_COUNT; i++) {
					cout << gBest[i] << ",";
				}
				average = 0;
				for (int i=0; i<SWARM_SIZE; i++) {
					average += evaluate(particles[i].x);
				}
				cout << "], AVERAGE: " << average/(float)iteration << endl;
				printCount *= 10;
			}

			//cout << "ITERATION: " << iteration << ", GLOBAL BEST: " << gBestResult << endl;
		}
	}

	return 0;
}



