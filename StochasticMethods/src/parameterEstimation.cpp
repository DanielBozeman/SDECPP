#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "RandomUtils.hpp"
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <limits>

struct KeyHasher{
    static int GetHashCodeForBytes(const char * bytes, int numBytes)
        {
        unsigned long h = 0, g;
        for (int i=0; i<numBytes; i++)
        {
            h = ( h << 4 ) + bytes[i];
            if (g = h & 0xF0000000L) {h ^= g >> 24;}
            h &= ~g;
        }
        return h;
        }
        static int GetHashForDouble(double v)
{
   return GetHashCodeForBytes((const char *)&v, sizeof(v));
}

std::size_t operator()(std::vector<double> const& vec) const
{
   int ret = 0;
   for (int i=0; i<vec.size(); i++) ret += ((i+1)*(GetHashForDouble(vec[i])));
   return ret;
}
};

double rmse(std::vector<double> actual, std::vector<double> prediction){

    auto squareError = [](double a, double b) {
        auto e = a-b;
        return e*e;
    };
    double sum = 0;

    try{
        for(int i = 0; i < actual.size(); i++){
            sum += ((actual[i] - prediction[i])*(actual[i] - prediction[i]));
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

    std::unordered_map<std::vector<double>, double, KeyHasher> rmsMap;

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;

    double bestRMS = std::numeric_limits<double>::infinity();

    while (temperature > temperatureLimit){
        for(int i = 0; i < stepsAtTemp; i++){
            if (rmsMap.find(currentModel.parameters) != rmsMap.end())
            {
                RMS = rmsMap[currentModel.parameters];
            }else{
                curApproximation = averageEulerMaruyama(currentModel, numSimulations);
            
                RMS = rmse(observations, curApproximation);
            }
            
            prob = acceptanceProbability(RMS, oldRMS, temperature);


            if (prob > randomGenerator.d01()){
                oldRMS = RMS;
                std::cout << "\nNew RMS: " << RMS;
            }

            if(RMS < bestRMS){
                bestModel = currentModel;
                bestRMS = RMS;
            }

            parameterNeighbor(currentModel.parameters, currentModel.parameterLimits, currentModel.parameterSteps);       
        }
        temperature *= coolingRate;
        std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;
}


