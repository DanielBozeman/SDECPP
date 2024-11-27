#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "RandomUtils.hpp"
#include <algorithm>
#include <math.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <limits>

//Root mean square error between two vectors
//This can probably be changed to just MSE
double rmse(std::vector<double>& simulation, std::vector<double>& actual){
    double sum = 0;

    try{
        for(int i = 0; i < actual.size(); i++){
            sum += ((actual[i] - simulation[i])*(actual[i] - simulation[i]));
        }
        if(sum < 0){
            throw std::runtime_error("NOT A NUMBER REACHED!");
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    double rmse = (sum / actual.size()); 

    return sqrt(rmse);
}

//Root mean square error between the average of a bunch of simuations and a vector of observations
//Probably could also be swapped out with MSE
long double multiVectorRMSE(std::vector<std::vector<double>>& simulations, std::vector<double>& actual){

    auto squareError = [](double a, double b) {
        auto e = a-b;
        return e*e;
    };

    std::vector<double> average = {};

    average.reserve(simulations[0].size());

    double averageVal;

    for(int i = 0; i < simulations[0].size(); i++){
        averageVal = 0;
        for(int j = 0; j < simulations.size(); j++){
            averageVal += simulations[j][i];
        }
        average.push_back(averageVal);
    }

    for(int i = 0; i < average.size(); i++){
        average[i] /= simulations.size();
    }

    return rmse(average, actual);
}

//Finds the set of returns of a set of simulated paths
//Return is organized as the first coordinate being the time and second being the particular return returns[timeStep][return]
std::vector<std::vector<double>> findReturns(std::vector<std::vector<double>>& simulations){

    std::vector<std::vector<double>> simReturns = {};

    for(int i = 1; i < simulations[0].size(); i++){
        simReturns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            simReturns[i-1].push_back(simulations[j][i] - simulations[j][i-1]);
        }
    }

    return simReturns;
}

//Returns the returns of a set of data at each timestep, but this time in absolute value
//Can be useful for sorting and may have better statistical properties
std::vector<std::vector<double>> findAbsReturns(std::vector<std::vector<double>>& simulations){

    std::vector<std::vector<double>> simReturns = {};

    for(int i = 1; i < simulations[0].size(); i++){
        simReturns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            simReturns[i-1].push_back(abs(simulations[j][i] - simulations[j][i-1]));
        }
    }

    return simReturns;
}

//Calculates the sample mean of a set of data
double sampleMean(std::vector<double>& observations){
    double mean = 0;

    int skips = 0;

    for(int i = 0; i < observations.size(); i++){
        if(!isnan(observations[i])){
            mean += observations[i]; 
        }else{
            skips += 1;
        }
    }

    //std::cout << "\nMean: " << mean;

    mean /= (observations.size() - skips);

    return mean;
}

//Calculates the sample variance of a set of data
double sampleVariance(std::vector<double>& observations){

    //double mean = sampleMean(observations);
    double mean = 0;

    double variance = 0;

    int skips=0;

    for(int i = 0; i < observations.size(); i++){
        if(!isnan(observations[i])){
            variance += (observations[i] - mean)*(observations[i] - mean);
            //std::cout << "\nVar: " << variance;
        }else{
            skips += 1;
        }
    }

    //std::cout << "\nVar: " << variance;
    variance /= (observations.size() - 1.0 - skips);


    return variance;
}

//Calculates the pdf of a normal distribution with given mean and variance at a given datapoint
long double normalPDF(double observation, double mean, double variance){
    long double pdf = 1/sqrt(2 * 3.1415926 * variance);

    long double endPDF = pdf * exp(-1* (((observation - mean)*(observation - mean))/(2 * variance)));


    //std::cout << "\npdf part: " << pdf;

    //std::cout << "\ne part: " << exp(-1* (((observation - mean)*(observation - mean))/(2 * variance)));

    return endPDF;
}

//Calculates the cdf of a normal distribution with given mean and variance at a given datapoint
double normalCDF(double& observation, double& mean, double& variance){

    double value = 0;
    if(observation > mean){
        value = 2 * mean - observation;
    }else{
        value = observation;
    }

    double cdf = (value- mean)/(sqrt(variance * 2));

    cdf = 0.5 * (1 + erf(cdf));

    //std::cout << "\nCDF is " << cdf;

    return cdf;
}

//COST FUNCTION
//Compares the returns of the given data against a set of simulations
//Should get smaller as the simulated variance gets more accurate
//Still in progress
long double returnComparison(std::vector<std::vector<double>>& simulations, std::vector<double>& actual){

    std::vector<std::vector<double>> simReturns = {};

    for(int i = 1; i < simulations[0].size(); i++){
        simReturns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            simReturns[i-1].push_back(simulations[j][i] - simulations[j][i-1]);
        }
    }

    long double totalCost = 1;

    for(int i = 0; i < (actual.size()); i++){
        //std::cout << "\nReturn size: " << simReturns.size();

        double simMean = sampleMean(simReturns[i]);
        simMean = 0;
        double simVariance = sampleVariance(simReturns[i]);

        //std::cout << "\nVar: " << simVariance;
        //std::cout << "\nMean: " << simMean;
        //std::cout << "\nActual: " << actual[i];

        long double pdf = normalPDF(actual[i], simMean, simVariance);

        //long double cdf = normalCDF(actual[i], simMean, simVariance);

        //std::cout << "\nCDF: " << cdf;

        totalCost += (log(pdf));
        
        //std::cout << "\nTotal: " << totalCost << " Log: " << log(cdf);
        //std::cout << "\nMean: " << simMean << " Var: " << simVariance << " Obs: " << actual[i];
    }

    return (-1 * totalCost);
}

//Acceptance Probability for simulated annealing
double acceptanceProbability(double newState, double oldState, double temperature){
    //std::cout << "\nNew: " << newState << "    Old: " << oldState;
    if (newState <= oldState){
        return 1;
    }
    else{
        return(exp(-1 * (newState - oldState) / temperature));
    }
}

//Finds a random neighbor of a parameter
std::vector<double> parameterNeighbor(std::vector<double> currentParameters, std::vector<std::vector<double>>& parameterLimits, std::vector<double>& parameterStepSize){
    int choice = randomGenerator.next() % currentParameters.size();

    currentParameters[choice] += randomPathMaker::dW(parameterStepSize[choice]);

    //std::cout << "\nVal: " << currentParameters[choice];

    currentParameters[choice] = currentParameters[choice] < parameterLimits[choice][0] ? parameterLimits[choice][0] : currentParameters[choice];
    currentParameters[choice] = currentParameters[choice] > parameterLimits[choice][1] ? parameterLimits[choice][1] : currentParameters[choice];

    return currentParameters;
}

//Completely randomizes a parameter
std::vector<double> randomParam(std::vector<double> currentParameters, std::vector<std::vector<double>>& parameterLimits){
    int choice = randomGenerator.next() % currentParameters.size();

    int limitSize = abs(parameterLimits[choice][1] - parameterLimits[choice][0]);

    currentParameters[choice] = (randomGenerator.next() % limitSize) + randomGenerator.d01() + std::min({parameterLimits[choice][1], parameterLimits[choice][0]});

    return currentParameters;
}

//Performs simulated annealing by changing a given parameter, tries to minimize the costFunction
std::vector<double> simulatedAnnealingVolEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost){
    double temperature = startingTemperature;

    std::vector<std::vector<double>> curApproximation;
    double RMS;
    double prob;
    double oldRMS = std::numeric_limits<double>::infinity();

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;
    stochasticModel newModel = model;

    stochasticModel noVarModel = model;
    noVarModel.betaFunction = zeroFunction;

    std::vector<double> drift = eulerMaruyama(noVarModel);

    for(int i = 0; i < observations.size(); i++){
        observations[i] = observations[i] - drift[i];
    }

    std::vector<double> trueReturns = {};

    trueReturns.reserve(observations.size() - 1);

    for(int i = 1; i < observations.size(); i++){
        trueReturns.push_back((observations[i] - observations[i-1]));
    }

    double bestRMS = std::numeric_limits<double>::infinity();

    while (temperature > temperatureLimit){

        bool moved = false;

        for(int i = 0; i < stepsAtTemp; i++){

            curApproximation = multipleEulerMaruyama(newModel, numSimulations);

            for(int j = 0; j < curApproximation.size(); j++){
                for(int k = 0; k < curApproximation[k].size(); k++){
                    curApproximation[j][k] = curApproximation[j][k] - drift[k];
                }
            }
            
            RMS = cost(curApproximation, trueReturns);

            //std::cout << "\nCurrent Var: " << newModel.parameters[1][0] << "    " << newModel.parameters[1][1] << "    " << newModel.parameters[1][2];
            //std::cout << "    Current Cost: " << RMS;

            if((RMS == 0 || RMS == std::numeric_limits<double>::infinity() || std::isnan(RMS)) && (oldRMS == 0 || oldRMS == std::numeric_limits<double>::infinity() || std::isnan(RMS))){
                newModel.parameters[parameterSet] = randomParam(newModel.parameters[parameterSet], newModel.parameterLimits[parameterSet]);
                continue;
            }else{
                prob = acceptanceProbability(RMS, oldRMS, temperature);
            }

            if (prob > randomGenerator.d01()){
                oldRMS = RMS;
                currentModel = newModel;
                moved = true;
                //std::cout << "\nCurrent Var: " << newModel.parameters[1][0] << "  " << newModel.parameters[1][1];
                //std::cout << "    New Cost: " << RMS;
            }

            if(RMS < bestRMS){
                bestModel = currentModel;
                bestRMS = RMS;
            }

            newModel.parameters[parameterSet] = parameterNeighbor(currentModel.parameters[parameterSet], currentModel.parameterLimits[parameterSet], currentModel.parameterSteps[parameterSet]);       

        }
        if(!moved){
            std::cout << "\nNo movement that step, done here";
            break;
        }
        temperature *= coolingRate;
        std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;

    return model.parameters[parameterSet];
}

std::vector<double> simulatedAnnealingDriftEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost){
    double temperature = startingTemperature;

    std::vector<std::vector<double>> curApproximation;
    double RMS;
    double prob;
    double oldRMS = std::numeric_limits<double>::infinity();

    model.betaFunction = zeroFunction;

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;
    stochasticModel newModel = model;

    double bestRMS = std::numeric_limits<double>::infinity();

    while (temperature > temperatureLimit){
        for(int i = 0; i < stepsAtTemp; i++){

            curApproximation = multipleEulerMaruyama(newModel, numSimulations);
            
            RMS = cost(curApproximation, observations);
            //std::cout << "\nCurrent Var: " << newModel.parameters[0][0] << "   " << newModel.parameters[0][1];
            //std::cout << "    Current Cost: " << RMS;

            if((RMS == 0 || RMS == std::numeric_limits<double>::infinity()) && (oldRMS == 0 || oldRMS == std::numeric_limits<double>::infinity())){
                newModel.parameters[parameterSet] = randomParam(currentModel.parameters[parameterSet], currentModel.parameterLimits[parameterSet]);
                continue;
            }else{
                prob = acceptanceProbability(RMS, oldRMS, temperature);
            }

            if (prob > randomGenerator.d01()){
                oldRMS = RMS;
                currentModel = newModel;
                //std::cout << "\nCurrent Var: " << newModel.parameters[1][0];
                //std::cout << "    New Cost: " << RMS;
            }

            if(RMS < bestRMS){
                bestModel = currentModel;
                bestRMS = RMS;
            }

            newModel.parameters[parameterSet] = parameterNeighbor(currentModel.parameters[parameterSet], currentModel.parameterLimits[parameterSet], currentModel.parameterSteps[parameterSet]);       
        }
        temperature *= coolingRate;
        std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;

    return model.parameters[parameterSet];
}


