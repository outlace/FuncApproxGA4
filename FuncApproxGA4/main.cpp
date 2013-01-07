//
//  main.cpp
//  FuncApproxGA4
//
//  Created by Brandon Brown on 7/29/11.
//  Copyright 2011 UCLA. All rights reserved.
//
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;

// Each term consists of a multiplier (coefficient) and an exponent (i.e. 3 terms = 3 multipliers and 3 exponents)
// Exponents and multipliers are separated in different arrays and only recombine with each other (i.e. exponent numbers only recombine
// with other exponent numbers.

#define initPopNum 500
#define numTerms 8
#define numGenerations 20000
#define mutationRate 0.005
class Solution;


void initPop();
void showInd(Solution);
double runFunction(Solution,double);
void evaluate(Solution*);
void normalizeScores();
void mate();
void run();
int ranNum(int);
static double curAvgFit;
bool compare(Solution, Solution);
double ranDoub();
void writeData();
///
static int genNum;
static int numSols;
class Solution {
public:
    int id;
    double exp[numTerms];
    double cos[numTerms];
    long double fitScore;
    Solution() {
        id = numSols++;
    }
};
// used to store avgFit scores each generation
double computedData[numGenerations][2]; // stores avg fitness as each generation passes. Use to see if algorithm is smoothly progressing
/// current generation of solutions:
Solution solutions[initPopNum];
/// sample data used for testing only:
double data[8][2];
/// MAIN ALGORITHM FUNCTIONS: ----

void initPop() {
    // Initialize random population of size initPopNum
    curAvgFit = 1.0;
    for (int i = 0; i < initPopNum; i++) {
        //Solution* a = new Solution();
        solutions[i] = *(new Solution());
        // Populate arrays with random values
        for (int z = 0; z < numTerms; z++) {
            solutions[i].exp[z] = (ranNum(10));
            solutions[i].cos[z] = ((double)ranNum(10));
        }
        //cout << solutions[i].id <<endl;
        //delete a;
    }
    
    
}

void showInd(Solution sol) {
    
    for (int i = 0; i < numTerms; i++ ) {
        cout << sol.cos[i] << "x^" << sol.exp[i];
        if (i != (numTerms - 1))
            cout << " + ";
    }
    
}

double runFunction(Solution sol, double x) {
    double y = 0;
    for (int i = 0; i < numTerms; i++) {
        y += (sol.cos[i] * pow(x,sol.exp[i]));
    }
    
    return y;
}
// evaluates a single solution and assign it a fitScore
void evaluate(Solution* sol) {
    long double leastSquares = 0;
    // iterate through all sample X values
    for (int i = 0; i < 7; i++) {
        double xval = data[i][0];
        double obsY = runFunction(*sol, xval);
        double expY = data[i][1];
        double residual = (obsY - expY);
        // square residual
        leastSquares += powf(residual, 2);
    }
    sol->fitScore = leastSquares;
}

void normalizeScores() {
    // Run evaluate on each solution
    double totalFit = 0;
    for (int i = 0; i < initPopNum; i++) {
        evaluate(&(solutions[i]));
        if (solutions[i].fitScore < 2) {
            cout << "Optimal solution found! " << "Fit Score: " << solutions[i].fitScore << " : "; showInd(solutions[i]);
            cout << endl;
            exit(0);
            break;
        }
        totalFit += solutions[i].fitScore;
    }
    double avgFit = (totalFit / (double)initPopNum);
    curAvgFit = avgFit;
    cout << "Avg Fitness (leastSquares) for generation: " << genNum << " is " << avgFit << endl;
    computedData[genNum][0] = genNum;
    computedData[genNum][1] = avgFit;
    // After evaluate() has been run on each solution, then the solutions will be sorted with lowest fitScores (leastSquares at top)
    // then the fitScores will be re-assigned based on their ranking
    //sort(solutions,solutions + sizeof(solutions)/sizeof(Solution), compare);
    // re-assign normalized (rank-based) fitScores
    for (int i = 0; i < initPopNum; i++) {
       // solutions[i].fitScore = (1 - ((double)i / (double)initPopNum));
        solutions[i].fitScore = log10l((double)solutions[i].fitScore + (double)1.0);
    }
    
}

void mate() {
    // next generation (temp storage)
    Solution nexGen[initPopNum];
    // tournament style selection process
    vector<int> ballsInBag;
    for (int g = 0; g < initPopNum; g++){
        /**/
        long double som = log2l(curAvgFit);
        long int sizeT2 = (double)(solutions[g].fitScore);
        long int sizeT = som*(1/(sizeT2+1.0));
        // Number of balls that should be added to the bag for this individual solution
        /*double thisFitScore = (100*solutions[g].fitScore);
        if (thisFitScore < 1) {
            thisFitScore = 1;
        }
        double sizeT = (curAvgFit * log10(thisFitScore));
        */
        if (sizeT < 1) {
            continue;
        }
        for (int h = 0; h < sizeT;h++) {
            ballsInBag.push_back(g);
            
        }
    }
    // select mating pairs then generate progeny
    for (int i = 0; i < (initPopNum / 2); i++) {
        int mate1 = ballsInBag[ranNum((int)ballsInBag.size())];
        int mate2 = ballsInBag[ranNum((int)ballsInBag.size())];
        // make sure mate1 != mate2
        if (mate1 == mate2)
            mate2 = ballsInBag[ranNum((int)ballsInBag.size())];
        
        int cutPt = ranNum(numTerms+1);
        if (cutPt == 0) {
            cutPt++;
        } else if (cutPt == (numTerms+1)) {
            cutPt--;
        }
        double exp_new1[numTerms], exp_new2[numTerms], cos_new1[numTerms], cos_new2[numTerms];
        // populate pieces before cutPt
        //cout << "Mate 1: "; showInd(solutions[mate1]); cout << endl;
        //cout << "Mate 2: "; showInd(solutions[mate2]); cout << endl;
        for (int a = 0; a < cutPt; a++) {
            exp_new1[a] = solutions[mate2].exp[a];
            cos_new1[a] = solutions[mate2].cos[a];
            exp_new2[a] = solutions[mate1].exp[a];
            cos_new2[a] = solutions[mate1].cos[a];
        }
        // finish populating after cutPt
        for (int c = cutPt; c < numTerms; c++) {
            exp_new1[c] = solutions[mate1].exp[c];
            cos_new1[c] = solutions[mate1].cos[c];
            exp_new2[c] = solutions[mate2].exp[c];
            cos_new2[c] = solutions[mate2].cos[c];
        }
       
        // populate nexGen
        for (int h = 0; h < numTerms; h++) {
            nexGen[i*2].exp[h] = exp_new1[h];
            nexGen[i*2].cos[h] = cos_new1[h];
            nexGen[(i*2)+1].exp[h] = exp_new2[h];
            nexGen[(i*2)+1].cos[h] = cos_new2[h];
        } 
        //cout << "Cut Point: " << cutPt << endl;
        //cout << "Mate_new1: "; showInd(nexGen[i*2]); cout << endl;
        //cout << "Mate_new2: "; showInd(nexGen[(i*2) + 1]); cout << endl << endl;
        
    }
        // Once all progeny have been produced, MUTATE
    for (int i = 0; i < initPopNum; i++) {
        if (ranDoub() < mutationRate) {
            // mutate by adding or multiplying by random number to a random value in either exp[] or cos[]
            int ranVal = ranNum(numTerms);
            
            if (ranDoub() >= 0.5) {
                // mutate exp[]:
                int mutNum = ranNum(2);
                if (ranDoub() >= 0.5) {
                    nexGen[i].exp[ranVal] = (nexGen[i].exp[ranVal] + mutNum);
                } else {
                    nexGen[i].exp[ranVal] = (nexGen[i].exp[ranVal] - mutNum);
                }
            } else {
                // else mutate cos[]:
                double mutNum = (double)(ranNum(2));
                    // add ranNum
                if (ranDoub() >= 0.5) {
                    nexGen[i].cos[ranVal] = (nexGen[i].cos[ranVal] + mutNum);
                } else {
                    nexGen[i].cos[ranVal] = (nexGen[i].cos[ranVal] - mutNum);
                }
            }
        }
    }
    
    // At this point progeny have been produced by recombination and mutated. Now to copy nexGen to current generation (solutions[])
    for (int r = 0; r < initPopNum; r++) {
        for (int f = 0; f < numTerms; f++) {
            solutions[r].exp[f] = nexGen[r].exp[f];
            solutions[r].cos[f] = nexGen[r].cos[f];
        }
    }
}


int ranNum (int max) {
    if (max > 0) {
        return (rand() % max);
    } else {
        return 0;
    }
}
// returns a double between 0 and 1
double ranDoub() {
    double ranDub = (((double)ranNum(10000000) / (double)10000000));
    return ranDub;
}

void run() {
        initPop();
    for (int i = 0; i < numGenerations; i++) {
        normalizeScores();
        mate();
        genNum++;
    }
    writeData();
}

bool compare(Solution sol1, Solution sol2) {
	return (sol1.fitScore < sol2.fitScore);
}

void writeData() {
    ofstream myfile;
    stringstream fileDest;
    fileDest << "/Users/brandonb/Desktop/FuncApproxGA_Data/FuncApproxGA-Data_" << ranNum(10000) << ".csv";
    string newString = fileDest.str();
    const char * blah = newString.c_str();
    myfile.open(blah);
    myfile << "Mutation Rate: " << mutationRate << " Num Generations: " << numGenerations << " Init Pop: " << initPopNum << " Num Terms: " << numTerms << endl;
    myfile << "GenNum, AvgFit" << endl;
    for (int ee = 0; ee < numGenerations; ee++) {
        myfile << computedData[ee][0] << "," << computedData[ee][1];
        myfile << endl;
    }
    myfile.close();
}

int main (int argc, const char * argv[])
{
    // make seed based on system time for unpredictable random numbers
    time_t seconds;
    time(&seconds);
    srand((unsigned int) seconds);
    // initialize numSols and genNum to 0
    genNum = 0;
    numSols = 0;
    // initialize sample data to model function: 4x^3 + 2x^2 + 3
    for (int i = 1; i < 20;i++) {
        for (int z = 0;z < 2;z++) {
            if (z == 0) {
                // set x value (0,1,2,3..)
                data[i-1][z] = i;
            } else {
                // set y value
                data[i-1][z] = (4*pow(i,3) + 2*pow(i,2) + 3);
            }
        }
    }
    // insert code here...
    cout << "Hello, World!\n";
    //Solution sol1, sol2;
    //cout << sol2.id << endl;;
    run();
    
    return 0;
}

