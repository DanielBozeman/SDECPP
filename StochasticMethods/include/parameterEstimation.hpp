#ifndef __PARAMETERESTIMATION_H__
#define __PARAMETERESTIMATION_H__

#include "RandomUtils.hpp"
#include "stochasticMethods.hpp"

typedef double (*costFunction)(std::vector<double>, std::vector<double>);

double rmse(std::vector<double> actual, std::vector<double> prediction);
double acceptanceProbability(double newState, double oldState, double temperature);
std::vector<double> parameterNeighbor(std::vector<double> currentParameters, std::vector<std::vector<double>> parameterLimits, std::vector<double> parameterStepSize);
std::vector<double> simulatedAnnealingParameterEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost);
#endif // __PARAMETERESTIMATION_H__