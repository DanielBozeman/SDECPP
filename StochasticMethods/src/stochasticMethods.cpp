#include "stochasticMethods.hpp"
#include "parameterEstimation.hpp"
#include <vector>
#include "RandomUtils.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

double zeroFunction(double& value, double& time, std::vector<double>& parameters){
    return 0;
}

double polynomialFunction(double& value, double& time, std::vector<double>& parameters){
    double total = 0;
    for(int i = 0; i < parameters.size(); i++){
        if(parameters[i] == 0){
            continue;
        }
        total += (parameters[i] * std::pow(value, i));
        double temp = std::pow(value, i);
        bool isNan = std::isnan(total);
        int temp2 = 0;
    }
    return total;
}

stochasticModel::stochasticModel(stochastic_function function1, stochastic_function function2, double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<std::vector<double>>> constantLimits, std::vector<std::vector<double>> stepSizes){
    alphaFunction = function1;
    betaFunction = function2;
    initialValue = startValue;
    timeInterval = times;
    parameters = constants;
    parameterLimits = constantLimits;
    parameterSteps = stepSizes;
}

void stochasticModel::setParameters(std::vector<std::vector<double>> constants){
    parameters = constants;
}

void stochasticModel::parameterNeighbor(int paramSet){
    int choice = randomGenerator.next() % parameters[paramSet].size();

    parameters[paramSet][choice] += randomPathMaker::dW(parameterSteps[paramSet][choice]);

    parameters[paramSet][choice] = parameters[paramSet][choice] < parameterLimits[paramSet][choice][0] ? parameterLimits[paramSet][choice][0] : parameters[paramSet][choice];
    parameters[paramSet][choice] = parameters[paramSet][choice] > parameterLimits[paramSet][choice][1] ? parameterLimits[paramSet][choice][1] : parameters[paramSet][choice];
}

void stochasticModel::randomizeParameter(int paramSet){
    int choice = randomGenerator.next() % parameters[paramSet].size();

    int limitSize = abs(parameterLimits[paramSet][choice][1] - parameterLimits[paramSet][choice][0]);

    parameters[paramSet][choice] = (randomGenerator.next() % limitSize) + randomGenerator.d01() + std::min({parameterLimits[paramSet][choice][1], parameterLimits[paramSet][choice][0]});
}

polynomialModel::polynomialModel(double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<int>> usedTerms, std::vector<std::vector<std::vector<double>>> constantLimits, std::vector<std::vector<double>> stepSizes):stochasticModel(polynomialFunction, polynomialFunction, startValue, times, constants, constantLimits, stepSizes){
    activeTerms = usedTerms;
    //stochasticModel(polynomialFunction, polynomialFunction, startValue, times, constants, constantLimits, stepSizes); 
}

void polynomialModel::parameterNeighbor(int paramSet){
    int choice = randomGenerator.next() % activeTerms[paramSet].size();

    parameters[paramSet][activeTerms[paramSet][choice]] += randomPathMaker::dW(parameterSteps[paramSet][choice]);

    parameters[paramSet][activeTerms[paramSet][choice]] = parameters[paramSet][activeTerms[paramSet][choice]] < parameterLimits[paramSet][activeTerms[paramSet][choice]][0] ? parameterLimits[paramSet][activeTerms[paramSet][choice]][0] : parameters[paramSet][activeTerms[paramSet][choice]];
    parameters[paramSet][activeTerms[paramSet][choice]] = parameters[paramSet][activeTerms[paramSet][choice]] > parameterLimits[paramSet][activeTerms[paramSet][choice]][1] ? parameterLimits[paramSet][activeTerms[paramSet][choice]][1] : parameters[paramSet][activeTerms[paramSet][choice]];
}

void polynomialModel::randomizeParameter(int paramSet){
    int choice = randomGenerator.next() % activeTerms[paramSet].size();

    int limitSize = abs(parameterLimits[paramSet][activeTerms[paramSet][choice]][1] - parameterLimits[paramSet][activeTerms[paramSet][choice]][0]);

    parameters[paramSet][activeTerms[paramSet][choice]] = (randomGenerator.next() % limitSize) + randomGenerator.d01() + std::min({parameterLimits[paramSet][activeTerms[paramSet][choice]][1], parameterLimits[paramSet][activeTerms[paramSet][choice]][0]});
}

void polynomialModel::addTerm(int paramSet, int term){
    if(std::find(activeTerms[paramSet].begin(), activeTerms[paramSet].end(), term) != activeTerms[paramSet].end()){
        return;
    }else{
        activeTerms[paramSet].push_back(term);
    }

    if(paramSet == 0){
        for(int i = parameters[paramSet].size(); i <= term; i++){
            parameterLimits[paramSet].push_back({-100,100});
            parameterSteps[paramSet].push_back(1);
            parameters[paramSet].push_back(0);
        }
    }else{
        parameters[paramSet] = {};
        parameterLimits[paramSet] = {};
        parameterSteps[paramSet] = {};
        for(int i = 0; i <= term; i++){
            parameterLimits[paramSet].push_back({0,100});
            parameterSteps[paramSet].push_back(1);
            parameters[paramSet].push_back(0);
        }
        for(int i = 0; i < activeTerms[paramSet].size(); i++){
            parameters[paramSet][activeTerms[paramSet][i]] = 0.01;
        }
    }

}

void polynomialModel::removeTerm(int paramSet, int term){
    auto place = std::find(activeTerms[paramSet].begin(), activeTerms[paramSet].end(), term);
    if(place == activeTerms[paramSet].end()){
        return;
        std::cout << "\nTerm not found!";
    }else{
        activeTerms[paramSet].erase(place);
        parameters[paramSet][term] = 0;
    }
}

void polynomialModel::removeLastTerm(int paramSet){
    int term = activeTerms[paramSet].back();
    removeTerm(paramSet, term);
}

void polynomialModel::removeAllTerms(int paramSet){
    for(int i = 0; i < activeTerms[paramSet].size(); i++){
        removeTerm(paramSet, activeTerms[paramSet][i]);
    }
}

void polynomialModel::removeAllTerms(){
    removeAllTerms(0);
    removeAllTerms(1);
}

void polynomialModel::addNextTerm(int paramSet){
    int term = *std::max_element(activeTerms[paramSet].begin(), activeTerms[paramSet].end());

    addTerm(paramSet, term+1);
}

void polynomialModel::addMultipleTerms(int paramSet, std::vector<int> terms){
    for(int i = 0; i < terms.size(); i++){
        addTerm(paramSet, terms[i]);
    }
}

void polynomialModel::addMultipleTerms(int paramSet, int maxTerm){
    for(int i = 0; i <= maxTerm; i++){
        addTerm(paramSet, i);
    }
}

void polynomialModel::setTermParameter(int paramSet, int term, double coefficient){
    parameters[paramSet][term] = coefficient;
}

std::vector<long double> polynomialModel::calculateAIC(std::vector<double>& observations, int numSims, int divisions, double percentage){

    std::vector<long double> likelihood = findLikelihood(*this, observations, numSims, divisions, percentage);

    return {((2 * (activeTerms[0].size() + activeTerms[1].size())) - 2*likelihood[0]), likelihood[1]};
}

void linearlySpacedVector(std::vector<double> &xs, double a, double b, double h){   
    std::size_t N = (b-a)/h + 1;
    xs.resize(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
} 

void linearlySpacedVectorBySize(std::vector<double> &xs, double a, double b, std::size_t N){
    double h = (b - a) / static_cast<double>(N);
    xs.resize(N+1);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
} 

std::vector<double> eulerMaruyama(stochasticModel model, std::vector<double> brownianPath)
{
    std::vector<double> approximation = {model.initialValue};

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPath.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPath = rp.makePath(model.timeInterval[0], model.timeInterval.back() + dt, dt);
    }

    approximation.reserve((model.timeInterval.back() - model.timeInterval[0])/dt);

    //std::cout << "\nPath size: " << brownianPath.size(); 

    for(int i = 1; i < model.timeInterval.size(); i++){
        double prevValue = approximation[i-1];
        double prevTime = model.timeInterval[i-1];
        double dW = brownianPath[i] - brownianPath[i-1];
        approximation.push_back((prevValue + model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt + model.betaFunction(prevValue, prevTime, model.parameters[1])*dW ));
    }

    return approximation;
}

std::vector<double> eulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, std::vector<double> brownianPath)
{
    std::vector<double> approximation = {initialValue};

    double dt = timeInterval[1] - timeInterval[0];

    if(brownianPath.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPath = rp.makePath(timeInterval[0], timeInterval.back(), dt);
    }

    approximation.reserve((timeInterval.back() - timeInterval[0])/dt);

    for(int i = 1; i < timeInterval.size(); i++){
        double prevValue = approximation[i-1];
        double prevTime = timeInterval[i-1];
        double dW = brownianPath[i] - brownianPath[i-1];
        approximation.push_back((prevValue + alphaFunction(prevValue, prevTime, parameters)*dt + betaFunction(prevValue, prevTime, parameters)*dW ));
    }

    return approximation;
}

void eulerMaruyamaByReference(std::vector<double> &approximation, stochasticModel model, std::vector<double> brownianPath){
    
    approximation.resize(model.timeInterval.size());

    approximation[0] = model.initialValue;

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPath.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPath = rp.makePath(model.timeInterval[0], model.timeInterval.back() + dt, dt);
    }

    for(int i = 1; i < model.timeInterval.size(); i++){
        double prevValue = approximation[i-1];
        bool prevNan = std::isnan(prevValue);
        double prevTime = model.timeInterval[i-1];
        double dW = brownianPath[i] - brownianPath[i-1];

        double nextDrift = prevValue + model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt;
        if(std::isinf(nextDrift)){
            approximation[i] = nextDrift;
            continue;
        }
        else{
            approximation[i] = nextDrift + model.betaFunction(prevValue, prevTime, model.parameters[1])*dW;
        }
        double nextAlpha = model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt;
        double nextBeta = model.betaFunction(prevValue, prevTime, model.parameters[1])*dW;
        double nextVal = (prevValue + model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt + model.betaFunction(prevValue, prevTime, model.parameters[1])*dW );
        bool nextNan = std::isnan(approximation[i]);
        //approximation[i] = ((prevValue + model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt + model.betaFunction(prevValue, prevTime, model.parameters[1])*dW ));
    }

    return;
}

void eulerMaruyamaWithin(std::vector<double> &approximation, stochasticModel model, int numDivisions, std::vector<double> brownianPath){
    std::vector<double> newTimes;

    linearlySpacedVectorBySize(newTimes, model.timeInterval[0], model.timeInterval.back(), (model.timeInterval.size()-1)*numDivisions);

    stochasticModel newModel = model;
    newModel.timeInterval = newTimes;

    std::vector<double> longApproximation;

    eulerMaruyamaByReference(longApproximation, newModel);

    approximation.resize(model.timeInterval.size());

    for(int i = 0; i < model.timeInterval.size()*numDivisions; i += numDivisions){
        approximation[i/numDivisions] = longApproximation[i];
        //std::cout << "\n" << newModel.timeInterval[i];
    }
}

void multipleEulerMaruyamaByReference(std::vector<std::vector<double>> &approximations, stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths){
    
    approximations.resize(numSimulations, std::vector<double>(model.timeInterval.size()));

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPaths.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPaths = rp.makeMultiplePaths(model.timeInterval[0], model.timeInterval.back() + dt, dt, numSimulations);
    }
   
    std::vector<double> approximation;

    for(int i = 0; i < numSimulations; i++){
        eulerMaruyamaByReference(approximations[i], model, brownianPaths[i]);
    }    
}

void multipleEulerMaruyamaWithin(std::vector<std::vector<double>> &approximations, stochasticModel model, int numDivisions, int numSimulations, std::vector<std::vector<double>> brownianPaths){

    approximations.resize(numSimulations, std::vector<double>(model.timeInterval.size()));

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPaths.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPaths = rp.makeMultiplePaths(model.timeInterval[0], model.timeInterval.back() + dt, dt, numSimulations);
    }
   
    std::vector<double> approximation;

    for(int i = 0; i < numSimulations; i++){
        
        eulerMaruyamaWithin(approximation, model, numDivisions, brownianPaths[i]);

        approximations[i] = approximation;
    }
}

std::vector<std::vector<double>> multipleEulerMaruyama(stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths){
    double prevValue;
    double prevTime;
    double dW;
    
    std::vector<std::vector<double>> approximations;

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPaths.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPaths = rp.makeMultiplePaths(model.timeInterval[0], model.timeInterval.back() + dt, dt, numSimulations);
    }

    approximations.reserve(((model.timeInterval.back() - model.timeInterval[0])/dt) * brownianPaths.size());

    for(int i = 0; i < numSimulations; i++){
        std::vector<double> approximation = {model.initialValue};

        approximation.reserve((model.timeInterval.back() - model.timeInterval[0])/dt);

        for(int j = 1; j < model.timeInterval.size(); j++){
            prevValue = approximation[j-1];
            prevTime = model.timeInterval[j-1];
            dW = brownianPaths[i][j] - brownianPaths[i][j-1];
            approximation.push_back((prevValue + model.alphaFunction(prevValue, prevTime, model.parameters[0])*dt + model.betaFunction(prevValue, prevTime, model.parameters[1])*dW ));
        }

        approximations.push_back(approximation);
    }

    return approximations;
}

std::vector<std::vector<double>> multipleEulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, int numSimulations, std::vector<std::vector<double>> brownianPaths){
    std::vector<std::vector<double>> approximations;

    double dt = timeInterval[1] - timeInterval[0];

    if(brownianPaths.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPaths = rp.makeMultiplePaths(timeInterval[0], timeInterval.back(), dt, numSimulations);
    }

    approximations.reserve(((timeInterval.back() - timeInterval[0])/dt) * brownianPaths.size());

    for(int i = 0; i < numSimulations; i++){
        std::vector<double> approximation = {initialValue};

        approximation.reserve((timeInterval.back() - timeInterval[0])/dt);

        for(int j = 1; j < timeInterval.size(); j++){
            double prevValue = approximation[j-1];
            double prevTime = timeInterval[j-1];
            double dW = brownianPaths[i][j] - brownianPaths[i][j-1];
            approximation.push_back((prevValue + alphaFunction(prevValue, prevTime, parameters)*dt + betaFunction(prevValue, prevTime, parameters)*dW ));
        }

        approximations.push_back(approximation);
    }

    return approximations;
}

std::vector<double> averageEulerMaruyama(stochasticModel model, int numSimulations){
    std::vector<std::vector<double>> approximations = multipleEulerMaruyama(model, numSimulations);

    std::vector<double> average = {};

    average.reserve(approximations[0].size());

    double averageVal;

    for(int i = 0; i < approximations[0].size(); i++){
        averageVal = 0;
        for(int j = 0; j < approximations.size(); j++){
            averageVal += approximations[j][i];
        }
        average.push_back(averageVal);
    }

    for(int i = 0; i < average.size(); i++){
        average[i] /= approximations.size();
    }

    return average;
}

double testFit(stochasticModel model, std::vector<double> data, int numInterRuns){
    
    stochasticModel baseModel = model;

    std::vector<double> times = model.timeInterval;

    double dt = ((times[1] - times[0])/1000);

    std::vector<std::vector<double>> endTests = {{}};
    
    for(int i = 1; i < times.size(); i++){

        std::vector<double> interTimes = {};

        interTimes.push_back(times[i-1]);

        int numTimes = (times[i] - times[i-1])/dt;

        for(int j = 0; j < numTimes; j++){
            interTimes.push_back(interTimes.back() + dt);
        }

        stochasticModel subModel = model;
        subModel.timeInterval = interTimes;
        subModel.initialValue = data[i-1];

        std::vector<double> finalValues = {};

        std::vector<std::vector<double>> interData = multipleEulerMaruyama(subModel, numInterRuns);

        for(int j = 0; j < interData.size(); j++){
            finalValues.push_back(interData[j].back());
        }

        std::sort(finalValues.begin(), finalValues.end());

        endTests.push_back(finalValues);
    }

    std::vector<int> ranks = {-5};

    for(int i = 1; i < data.size(); i++){
        auto it = std::lower_bound(endTests[i].begin(), endTests[i].end(), data[i]);

        int position = 0;
        if (it == endTests[i].end()) {
            //std::cout << "\nElement not found\n";
            //std::cout << "\nThe number " << data[i] << " would be inserted at position " << endTests[i].size() << std::endl;
            position = endTests[i].size();
        } else {
            position = std::distance(endTests[i].begin(), it);
            //std::cout << "\nThe number " << data[i] << " would be inserted at position " << position << std::endl;
        }

        ranks.push_back(position);
    }

    std::vector<int> omegas = {};

    for(int i = 0; i < numInterRuns + 1; i++){
        omegas.push_back(count(ranks.begin(), ranks.end(), i));
    }

    double chiVar = (double(data.size() - 1)/(numInterRuns + 1));

    //std::cout << "\n" << chiVar;

    double chiStatistic = 0;

    for(int i = 0; i < omegas.size(); i++){
        double temp = ((omegas[i] - chiVar)*(omegas[i]-chiVar))/double(chiVar);

        //std::cout << "\nCur stat: " << temp;

        chiStatistic += temp;
    }

    //std::cout << "\nChi var is " << chiVar;
    return chiStatistic;
}  
