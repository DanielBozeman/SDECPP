#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "RandomUtils.hpp"
#include "FileUtils.hpp"
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

long double multiVectorMSE(std::vector<std::vector<double>>& simulations, std::vector<double>& actual){
    
    auto squareError = [](double a, double b){
        auto e = a-b;
        return e*e;
    };

    double averageVal;

    long double sum = 0;

    for(int i = 0; i < actual.size(); i++){
        averageVal = 0;
        for(int j = 0; j < simulations.size(); j++){
            averageVal += simulations[j][i];
        }

        sum += squareError((averageVal/simulations.size()), actual[i]);
    }

    return (sum/actual.size());
    
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









//--------------------------------------------------------------------------------------------



/** -----------------------------------------------------------
 * @brief Cost function used for the drift function, finds MSE between model and real values
 * 
 * @param model SDE model to be used in the drift estimation
 * @param observations Data to be fit against 
 * @param numSims Number of simulations to use in the test, higher is more accurate
 * @param params Optional parameters, not used
 * @return long double MSE of model against the real data 
--------------------------------------------------------------- */
long double driftCost(stochasticModel model, std::vector<double> &observations, int numSims, std::vector<double> params){

    stochasticModel testModel = model;
    testModel.betaFunction = zeroFunction;

    std::vector<std::vector<double>> approximations;

    multipleEulerMaruyamaWithin(approximations, testModel,  params[0], numSims);

    long double mse = multiVectorMSE(approximations, observations);

    bool isNan = std::isnan(mse);

    return mse;
}

long double varianceCost(stochasticModel model, std::vector<double>& observations, int numSims, std::vector<double> params){

   stochasticModel noVarModel = model;
   noVarModel.betaFunction = zeroFunction;
   
   std::vector<double> noVarData;
   eulerMaruyamaWithin(noVarData, noVarModel, params[0]);

   std::vector<double> realVariance = {};

   for(int i = 1; i < observations.size(); i++){
        realVariance.push_back((observations[i] - noVarData[i]) - (observations[i-1] - noVarData[i-1]));
   }

   std::vector<std::vector<double>> simData;
   std::vector<double> simVar;

   stochasticModel simModel = model;

   double chance = 0;

   for(int i = 1; i < observations.size(); i++){
        simVar = {};

        double start = model.timeInterval[i-1];
        double end = model.timeInterval[i];
        std::vector<double> interTimes;
        linearlySpacedVectorBySize(interTimes, start, end, params[0]);

        simModel.timeInterval = interTimes;
        simModel.initialValue = observations[i-1];

        multipleEulerMaruyamaByReference(simData, simModel, params[1]);

        for(int j = 0; j < simData.size(); j++){
            simVar.push_back((simData[j].back() - noVarData[i]) - (simData[j][0] - noVarData[i-1]));
        }

        double var = sampleVariance(simVar);
        double mean = sampleMean(simVar);

        //std::cout << "\nVar: " << var << "\nMean: " << mean;

        double pdf = normalPDF(realVariance[i-1], mean, var);

        //std::cout << "\nPDF: " << pdf;

        chance += log(pdf);
   }

   return -1 * chance;
}

/** -----------------------------------------------------------
 * @brief Estimates the parameters of the alpha or beta function of a SDE model  
 * 
 * @param model Stochastic model with unkown parameters
 * @param parameterSet Parameter set to be estimated, 0 for alpha parameters, 1 for beta parameters
 * @param observations Data to be fit against
 * @param numSimsPerStep Number of simulations for each evalution of the cost function 
 * @param startingTemp Starting temperature of the simulated annealing algorithm    
 * @param coolingRate Rate at which the temperature drops at each step of the algorithm
 * @param stepsAtTemp Iterations of random parameters checked at each temperature   
 * @param tempLimit Lowest acceptable temperature, algorithm will stop at this temperature  
 * @param costFunction Cost function to be minimized 
 * @param optionalParams Optional parameters for the cost function 
 * @return std::vector<double> Fit parameters for the selected data set
--------------------------------------------------------------- */
std::vector<double> paramEstimation(stochasticModel model, int parameterSet, std::vector<double> observations, int numSimsPerStep, double startingTemp, double coolingRate, int stepsAtTemp, double tempLimit, modelCostFunction costFunction, std::vector<double> optionalParams){
    
    int moveLimit = 500;

    int notMovedIn = 0;
    
    double temperature = startingTemp;

    double cost;
    double prob;
    double oldCost = std::numeric_limits<double>::infinity();

    stochasticModel currentModel = model;
    stochasticModel bestModel = model;
    stochasticModel newModel = model;

    double bestCost = std::numeric_limits<double>::infinity();

    while(temperature > tempLimit){
        for(int i = 0; i < stepsAtTemp; i++){

            notMovedIn++;

            cost = costFunction(newModel, observations, numSimsPerStep, optionalParams);

            //std::cout << "\nCur cost: " << cost;
            //std::cout << "\nCur param" << newModel.parameters[0][0];

            if((cost == 0 || abs(cost) == std::numeric_limits<double>::infinity()) && (oldCost == 0 || abs(oldCost) == std::numeric_limits<double>::infinity())){
                newModel.randomizeParameter(parameterSet);
                //std::cout << "\nIndeterminate cost!";
                continue;
            }else{
                prob = acceptanceProbability(cost, oldCost, temperature);  
                //std::cout << "\nCur prob: " << prob;     
            }

            if(prob > randomGenerator.d01()){
                oldCost = cost;
                currentModel = newModel;
            }

            if(cost < bestCost){
                notMovedIn = 0;
                //std::cout << "\nBest cost: " << cost;
                bestModel = currentModel;
                bestCost = cost;
            }

            if(notMovedIn > moveLimit){
                return bestModel.parameters[parameterSet];
            }

            newModel.parameterNeighbor(parameterSet);
            //currentModel.parameterNeighbor(parameterSet);
            //newModel.parameters[parameterSet] = currentModel.parameters[parameterSet];
        }

        temperature *= coolingRate;
        //std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;
    
    return model.parameters[parameterSet];
}

std::vector<double> polynomialParamEstimation(polynomialModel model, int parameterSet, std::vector<double> observations, int numSimsPerStep, double startingTemp, double coolingRate, int stepsAtTemp, double tempLimit, modelCostFunction costFunction, std::vector<double> optionalParams){
    
    int moveLimit = 1000;

    int notMovedIn = 0;
    
    double temperature = startingTemp;

    double cost;
    double prob;
    double oldCost = std::numeric_limits<double>::infinity();

    polynomialModel currentModel = model;
    polynomialModel bestModel = model;
    polynomialModel newModel = model;

    double bestCost = std::numeric_limits<double>::infinity();

    while(temperature > tempLimit){
        for(int i = 0; i < stepsAtTemp; i++){

            notMovedIn++;

            cost = costFunction(newModel, observations, numSimsPerStep, optionalParams);

            //std::cout << "\nCur cost: " << cost;
            //std::cout << "\nNot moved in: " << notMovedIn; 
            //std::cout << "\nCur param" << newModel.parameters[0][0];

            if((cost == 0 || abs(cost) == std::numeric_limits<double>::infinity()) && (oldCost == 0 || abs(oldCost) == std::numeric_limits<double>::infinity())){
                newModel.randomizeParameter(parameterSet);
                //std::cout << "\nIndeterminate cost!";
                continue;
            }else{
                prob = acceptanceProbability(cost, oldCost, temperature);  
                //std::cout << "\nCur prob: " << prob;     
            }

            if(prob > randomGenerator.d01()){
                oldCost = cost;
                currentModel = newModel;
            }

            if(cost < bestCost){
                notMovedIn = 0;
                //std::cout << "\nBest cost: " << cost;
                bestModel = currentModel;
                bestCost = cost;
            }

            if(notMovedIn > moveLimit){
                return bestModel.parameters[parameterSet];
            }

            newModel.parameterNeighbor(parameterSet);
            //currentModel.parameterNeighbor(parameterSet);
            //newModel.parameters[parameterSet] = currentModel.parameters[parameterSet];
        }

        temperature *= coolingRate;
        //std::cout << "\nTemperature: " << temperature;
    }

    model = bestModel;
    
    return model.parameters[parameterSet];
}

polynomialModel polynomialMultiEstimation(polynomialModel model, std::vector<int> parameterSets, std::vector<double> observations, int numSimsPerStep, double startingTemp, double coolingRate, int stepsAtTemp, double tempLimit, modelCostFunction costFunction, std::vector<double> optionalParams){

    int moveLimit = 1000;

    int notMovedIn = 0;
    
    double temperature = startingTemp;

    double cost;
    double prob;
    double oldCost = std::numeric_limits<double>::infinity();

    polynomialModel currentModel = model;
    polynomialModel bestModel = model;
    polynomialModel newModel = model;

    double bestCost = std::numeric_limits<double>::infinity();

    while(temperature > tempLimit){
        for(int i = 0; i < stepsAtTemp; i++){

            notMovedIn++;

            cost = costFunction(newModel, observations, numSimsPerStep, optionalParams);

            //std::cout << "\nCur cost: " << cost;
            //std::cout << "\nNot moved in: " << notMovedIn; 
            //std::cout << "\nCur param" << newModel.toString();

            if((cost == 0 || abs(cost) == std::numeric_limits<double>::infinity()) && (oldCost == 0 || abs(oldCost) == std::numeric_limits<double>::infinity())){
                int choice = randomGenerator.next() % parameterSets.size();
                newModel.randomizeParameter(parameterSets[choice]);
                //std::cout << "\nIndeterminate cost!";
                continue;
            }else{
                prob = acceptanceProbability(cost, oldCost, temperature);  
                //std::cout << "\nCur prob: " << prob;     
            }

            if(prob > randomGenerator.d01()){
                oldCost = cost;
                currentModel = newModel;
            }

            if(cost < bestCost){
                notMovedIn = 0;
                //std::cout << "\nBest cost: " << cost << "\nBest model: \n" << bestModel.toString();
                bestModel = currentModel;
                bestCost = cost;
            }

            if(notMovedIn > moveLimit){
                return bestModel;
            }

            int choice = randomGenerator.next() % parameterSets.size();
            newModel.parameterNeighbor(parameterSets[choice]);
            //newModel.randomizeParameter(parameterSets[choice]);
            //currentModel.parameterNeighbor(parameterSet);
            //newModel.parameters[parameterSet] = currentModel.parameters[parameterSet];
        }

        temperature -= coolingRate;
        //std::cout << "\nTemperature: " << temperature;
    }

    return bestModel; 

}


double dtByPercentage(std::vector<double>& observations, double percentage, double input){


    std::sort(observations.begin(), observations.end());

    double dt = (observations.back() - observations[0])/observations.size();

    int numInRange = std::floor(observations.size() * percentage);
    
    std::vector<double>::iterator low = std::lower_bound(observations.begin(), observations.end(), input);
    
    int index = std::distance(observations.begin(), low);

    double max = observations.back();
    double min = observations[0];

    if(index == observations.size() || index == 0){
        return 0;
    }

    int upperIndex = index + std::floor(percentage*0.5*observations.size());
    if(upperIndex > observations.size() - 1){
        upperIndex = observations.size() - 1;
    }

    int lowerIndex = index - std::floor(percentage*0.5*observations.size());
    if(lowerIndex < 0){
        lowerIndex = 0;
    }

    double upperDistance = observations[upperIndex] - observations[index];
    double lowerDistance = observations[index] - observations[lowerIndex];

    //std::cout << "\nUpper: " << upperDistance << "   Lower: " << lowerDistance;

    return std::max(upperDistance,lowerDistance);
}

double findDt(std::vector<double>& observations, int divisions){
    double max = -std::numeric_limits<double>::infinity();
    double min = std::numeric_limits<double>::infinity();

    for(int i = 0; i < observations.size(); i++){
        if(observations[i] < min){
            min = observations[i];
        }
        if(observations[i] > max){
            max = observations[i];
        }
    }

    //std::cout << "\nMax: " << max;
    //std::cout << "\nMin: " << min;
    double dt = (max - min)/double(divisions);

    return dt;
}

double estimatePdf(std::vector<double>& observations, double input, double percentage){
    //double dt = findDt(observations, divisions);

    double dt = dtByPercentage(observations, percentage, input);

    if(dt == 0){
        return 0;
    }

    //dt = 0.05;
    double simLeft = input - dt;
    double simRight = input + dt;
    
    int countLeft = 0;
    int countRight = 0;

    
    for(int j = 0; j < observations.size(); j++){
        if(observations[j] <= simRight){
            countRight++;
        }
        if(observations[j] <= simLeft){
            countLeft++;
        }
    }

    double dy = (double(countRight)/observations.size()) - (double(countLeft)/observations.size());
    double dx = simRight - simLeft;

    double pdf = dy/dx;

    //std::cout << "\nCount Left: " << countLeft;
    //std::cout << "\nCount Right: " << countRight;

    //std::cout << "\ndy: " << dy;
    //std::cout << "\ndx: " << dx;

    return pdf;
}

std::vector<long double> findLikelihood(stochasticModel model, std::vector<double> observations, int numSims, int divisions, double percentage){  

   std::vector<long double> returnVector = {0, 0};

   std::vector<std::vector<double>> simData;
   std::vector<double> simEnds;

   stochasticModel simModel = model;

   long double chance = 0;

   long int exponent = 0;
   long double mantissa = 0;

   int tempExp = 0;

   multipleEulerMaruyamaWithin(simData, simModel, divisions, numSims);

   for(int i = 1; i < observations.size(); i++){
        simEnds = {};

        for(int j = 0; j < simData.size(); j++){
            double data = simData[j][i];
            simEnds.push_back(simData[j][i]);
        }
        
        double pdf = estimatePdf(simEnds, observations[i], percentage);
    
        // double tempMantissa = std::frexp(pdf, &tempExp);

        // mantissa += log2(tempMantissa);
        // exponent += tempExp;

        // exponent += floor(mantissa);
        // mantissa -= floor(mantissa);

        returnVector[0] += log(pdf);

        if(std::isinf(returnVector[0])){
            returnVector[1] = i;
            return returnVector;
        }
   }

   //std::cout << "\nMantissa: " << mantissa << "   Exp:" << exponent;

   return returnVector;
}

//RETURNS TRUE IF AIC1 IS SMALLER THAN AIC2
bool compareAIC(std::vector<long double> AIC1, std::vector<long double> AIC2){

    //If both AICs are not inf then we just compare them
    if(!std::isinf(AIC1[0]) && !std::isinf(AIC2[0])){
        return(AIC1[0] < AIC2[0]);
    }

    //If both are inf then we compare how far they worked
    if(std::isinf(AIC1[0]) && std::isinf(AIC2[0])){
        return(AIC1[1] > AIC2[1]);
    }

    //If one AIC is inf but not the other then we can return true or false depending
    if(std::isinf(AIC1[0])){
        return false;
    }
    else if(std::isinf(AIC2[0])){
        return true;
    }

    return false;
}