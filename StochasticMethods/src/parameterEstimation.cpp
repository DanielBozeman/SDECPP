#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "RandomUtils.hpp"
#include <algorithm>
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <limits>

double rmse(std::vector<double> simulation, std::vector<double> actual){
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
    
    double rmse = sqrt((sum / actual.size())); 

    return rmse;
}

double multiVectorRMSE(std::vector<std::vector<double>> simulations, std::vector<double> actual){

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

double findMax(std::vector<double> dataVector){
    
    double maximum = 0;

    for(int i = 0; i < dataVector.size(); i++){
        if(abs(dataVector[i]) > maximum){
            maximum = abs(dataVector[i]);
        }
    }

    return maximum;
}

double averageLogReturnComparison(std::vector<std::vector<double>> simulations, std::vector<double> actual){

    std::vector<std::vector<double>> returns;

    for(int i = 1; i < simulations[0].size(); i++){
        returns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            returns[i-1].push_back(log(simulations[j][i]/simulations[j][i-1]));
        }
    }

    double logDifference = 0;

    std::vector<double> logReturns = {};

    logReturns.reserve(actual.size() - 1);

    for(int i = 1; i < actual.size(); i++){
        logReturns.push_back(log(actual[i]/actual[i-1]));\
    }
 
    double difference = 0;

    for(int i = 0; i < logReturns.size(); i++){
        difference = abs(logReturns[i] - findMax(returns[i]));

        logDifference += difference;
    }

    return logDifference;
}

double sampleMean(std::vector<double> observations){
    double mean = 0;

    for(int i = 0; i < observations.size(); i++){
        mean += observations[i];
    }

    mean /= observations.size();

    return mean;
}

double sampleVariance(std::vector<double> observations){

    double mean = sampleMean(observations);

    double variance = 0;

    for(int i = 0; i < observations.size(); i++){
        variance += (observations[i] - mean)*(observations[i] - mean);
    }

    variance /= (observations.size() - 1.0);

    return variance;
}

double normalPDF(double observation, double mean, double variance){
    double pdf = 1/sqrt(2 * 3.1415926 * variance);

    pdf *= exp(-1 * (((observation - mean)*(observation - mean))/(2 * variance)));

    return pdf;
}

double normalCDF(double observation, double mean, double variance){

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

double returnComparison(std::vector<std::vector<double>> simulations, std::vector<double> actual){
    std::vector<std::vector<double>> simReturns;

    for(int i = 1; i < simulations[0].size(); i++){
        simReturns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            simReturns[i-1].push_back(simulations[j][i] - simulations[j][i-1]);
        }
    }

    std::vector<double> trueReturns = {};

    trueReturns.reserve(actual.size() - 1);

    for(int i = 1; i < actual.size(); i++){
        trueReturns.push_back(actual[i] - actual[i-1]);
    }

    double totalCost = 0;

    double simMean = 0;
    double simVariance = 0;

    simMean = sampleMean(simReturns[10]);
    simVariance = sampleVariance(simReturns[10]);

    double pdf = normalPDF(trueReturns[10], simMean, simVariance);

    for(int i = 0; i < trueReturns.size(); i++){
         simMean = sampleMean(simReturns[i]);
         simVariance = sampleVariance(simReturns[i]);

         double pdf = normalPDF(trueReturns[i], simMean, simVariance);

         totalCost -= pdf;
    }

    return totalCost;
}

double acceptanceProbability(double newState, double oldState, double temperature){
    if (newState < oldState){
        return 1;
    }
    else{
        return(exp(-1 * (newState - oldState) / temperature));
    }
}

std::vector<double> parameterNeighbor(std::vector<double> currentParameters, std::vector<std::vector<double>> parameterLimits, std::vector<double> parameterStepSize){

    int upDown = (randomGenerator.next()%2)*2 - 1;

    int choice = randomGenerator.next() % currentParameters.size();

    //currentParameters[choice] += upDown * parameterStepSize[choice];

    currentParameters[choice] += randomPathMaker::dW(parameterStepSize[choice]);

    currentParameters[choice] = currentParameters[choice] < parameterLimits[choice][0] ? parameterLimits[choice][0] : currentParameters[choice];
    currentParameters[choice] = currentParameters[choice] > parameterLimits[choice][1] ? parameterLimits[choice][1] : currentParameters[choice];

    return currentParameters;
}

std::vector<double> simulatedAnnealingParameterEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit, costFunction cost){
    double temperature = startingTemperature;

    std::vector<std::vector<double>> curApproximation;
    double RMS;
    double prob;
    double oldRMS = std::numeric_limits<double>::infinity();

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;
    stochasticModel newModel = model;

    double bestRMS = std::numeric_limits<double>::infinity();

    while (temperature > temperatureLimit){
        for(int i = 0; i < stepsAtTemp; i++){

            curApproximation = multipleEulerMaruyama(newModel, numSimulations);
            
            RMS = cost(curApproximation, observations);
            //std::cout << "\nCurrent Cost: " << RMS;

            prob = acceptanceProbability(RMS, oldRMS, temperature);

            if (prob > randomGenerator.d01()){
                oldRMS = RMS;
                currentModel = newModel;
                //std::cout << "\nNew Cost: " << RMS;
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


