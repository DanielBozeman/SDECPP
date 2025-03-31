#ifndef __PARAMETERESTIMATION_H__
#define __PARAMETERESTIMATION_H__

#include "RandomUtils.hpp"
#include "stochasticMethods.hpp"

typedef long double (*costFunction)(std::vector<std::vector<double>>&, std::vector<double>&);

typedef long double (*modelCostFunction)(stochasticModel model, std::vector<double>&, int, std::vector<double>);



long double multiVectorRMSE(std::vector<std::vector<double>>& simulations, std::vector<double>& actual);
long double multiVectorMSE(std::vector<std::vector<double>>& simulations, std::vector<double>& actual);

double sampleMean(std::vector<double>& observations);
double sampleVariance(std::vector<double>& observations);
long double normalPDF(double observation, double mean, double variance);
double normalCDF(double& observation, double& mean, double& variance);
long double returnComparison(std::vector<std::vector<double>>& simulations, std::vector<double>& actual);
double averageLogReturnComparison(std::vector<std::vector<double>> simulations, std::vector<double> actual);
double acceptanceProbability(double newState, double oldState, double temperature);
std::vector<double> simulatedAnnealingVolEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost);
std::vector<double> simulatedAnnealingDriftEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost);

long double driftCost(stochasticModel model, std::vector<double> &observations, int numSims, std::vector<double> optionalParams = {});
long double varianceCost(stochasticModel model, std::vector<double> &observations, int numSims, std::vector<double> optionalParams = {});
std::vector<double> paramEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimsPerStep, double startingTemp, double coolingRate, int stepsAtTemp, double tempLimit, modelCostFunction costFunction, std::vector<double> optionalParams = {});

double findLikelihood(stochasticModel model, std::vector<double> observations, int numSims, int divisions, double dt);

#endif // __PARAMETERESTIMATION_H__