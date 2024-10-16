#ifndef BROWNIANPATHS
#define BROWNIANPATHS


#include <vector>
#include "Xoshiro256plus.h"

Xoshiro256plus timeSeedRand();

extern Xoshiro256plus randomGenerator;

class randomPathMaker{
    private:
    static double inv_normal_cdf(double);
    public:
    randomPathMaker();
    static double dW(double);
    std::vector<double> makePath(double intervalStart, double intervalEnd, double timeStep);
    std::vector<std::vector<double>> makeMultiplePaths(double intervalStart, double intervalEnd, double timeStep, int numPaths);
    std::vector<std::vector<double>> makeCorrelatedPaths(double, double, double, double);
};

#endif