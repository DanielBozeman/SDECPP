#ifndef __PARAMETERESTIMATION_H__
#define __PARAMETERESTIMATION_H__

#include "RandomUtils.hpp"
#include "stochasticMethods.hpp"

double rmse(std::vector<double> actual, std::vector<double> prediction);
double acceptanceProbability(double newState, double oldState, double temperature);
std::vector<double> parameterNeighbor(std::vector<double> currentParameters, std::vector<std::vector<double>> parameterLimits, std::vector<double> parameterStepSize);
void simulatedAnnealingParameterEstimation(stochasticModel& model, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit);
#endif // __PARAMETERESTIMATION_H__