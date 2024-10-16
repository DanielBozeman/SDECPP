#include "stochasticMethods.hpp"
#include <vector>
#include "RandomUtils.hpp"

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

std::vector<double> eulerMaruyama(stochasticModel model, std::vector<double> brownianPath)
{
    std::vector<double> approximation = {model.initialValue};

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPath.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPath = rp.makePath(model.timeInterval[0], model.timeInterval.back(), dt);
    }

    approximation.reserve((model.timeInterval.back() - model.timeInterval[0])/dt);

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

std::vector<std::vector<double>> multipleEulerMaruyama(stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths){
    double prevValue;
    double prevTime;
    double dW;
    
    std::vector<std::vector<double>> approximations;

    double dt = model.timeInterval[1] - model.timeInterval[0];

    if(brownianPaths.size() == 0){
        randomPathMaker rp = randomPathMaker();
        brownianPaths = rp.makeMultiplePaths(model.timeInterval[0], model.timeInterval.back(), dt, numSimulations);
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
