#ifndef __MODELESTIMATIONTESTS_H__
#define __MODELESTIMATIONTESTS_H__

    #include "stochasticMethods.hpp"
    #include <string>

    stochasticModel fitBlackScholes(std::string fileName, int dataColumn, int dataStart = 0, int dataEnd = 0);
    stochasticModel fitRandom(std::string fileName, int dataColumn, int dataStart = 0, int dataEnd = 0);
    stochasticModel fitOrnstein(std::string fileName, int dataColumn, int dataStart = 0, int dataEnd = 0);
    stochasticModel fitCEV();

    polynomialModel fitPolynomial(std::vector<double> &observations, std::vector<double> times, int maxTerm, int divisions);
    
    void createPath();
#endif // __MODELESTIMATIONTESTS_H__