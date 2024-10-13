#ifndef __PARAMETERESTIMATION_H__
#define __PARAMETERESTIMATION_H__

#include "RandomUtils.hpp"

double rmse(std::vector<double> actual, std::vector<double> prediction);
double acceptanceProbability(double newState, double oldState, double temperature);
void parameterNeighbor(std::vector<double>& currentParameters, std::vector<std::vector<double>> parameterLimits, std::vector<double> parameterStepSize);

#endif // __PARAMETERESTIMATION_H__