#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "RandomUtils.hpp"
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <limits>



double rmse(std::vector<double> actual, std::vector<double> prediction){

    auto squareError = [](double a, double b) {
        auto e = a-b;
        return e*e;
    };

    double sum = std::transform_reduce(actual.begin(), actual.end(), prediction.begin(), 0, std::plus<>(), squareError);
    double rmse = sqrt((sum / actual.size())); 

    return rmse;
}

double acceptanceProbability(double newState, double oldState, double temperature){
    if (newState < oldState){
        return 1;
    }
    else{
        return(exp(-1 * (newState - oldState) / temperature));
    }
}



void parameterNeighbor(std::vector<double>& currentParameters, std::vector<std::vector<double>> parameterLimits, std::vector<double> parameterStepSize){

    int upDown = (randomGenerator.next()%2)*2 - 1;

    int choice = randomGenerator.next() % currentParameters.size();

    currentParameters[choice] += upDown * parameterStepSize[choice];

    currentParameters[choice] = currentParameters[choice] < parameterLimits[choice][0] ? parameterLimits[choice][0] : currentParameters[choice];
    currentParameters[choice] = currentParameters[choice] > parameterLimits[choice][1] ? parameterLimits[choice][1] : currentParameters[choice];

    return;
}

void simulatedAnnealingParameterEstimation(stochasticModel& model, std::vector<double> observations, int numSimulations, double startingTemperature, double coolingRate, int stepsAtTemp, double temperatureLimit){
    double temperature = startingTemperature;

    std::vector<double> curApproximation;
    double RMS;
    double prob;
    double oldRMS = std::numeric_limits<double>::infinity();

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;

    while (temperature > temperatureLimit){
        for(int i = 0; i < stepsAtTemp; i++){
            curApproximation = averageEulerMaruyama(currentModel, numSimulations);

            
            RMS = rmse(observations, curApproximation);

            prob = acceptanceProbability(RMS, oldRMS, temperature);

            if (prob > randomGenerator.d01()){
                bestModel = currentModel;
            }

            parameterNeighbor(currentModel.parameters, currentModel.parameterLimits, currentModel.parameterSteps);       
        }
        temperature *= coolingRate;
        std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;
}


