#ifndef __STOCHASTICMETHODS_H__
#define __STOCHASTICMETHODS_H__

#include <vector>

typedef double (*stochastic_function)(double&, double&, std::vector<double>&);

class stochasticModel{
    private:

    public:
    stochastic_function alphaFunction;
    stochastic_function betaFunction;
    double initialValue;
    std::vector<double> timeInterval;
    std::vector<std::vector<double>> parameters;
    std::vector<std::vector<std::vector<double>>> parameterLimits;
    std::vector<std::vector<double>> parameterSteps;
    stochasticModel(stochastic_function function1, stochastic_function function2, double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<std::vector<double>>> constantLimits = {}, std::vector<std::vector<double>> stepSizes = {});  
    void setParameters(std::vector<std::vector<double>> constants);  
};

double zeroFunction(double& value, double& time, std::vector<double>& parameters);

std::vector<double> eulerMaruyama(stochasticModel model, std::vector<double> brownianPath = {});

std::vector<double> eulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, std::vector<double> brownianPath = {});

std::vector<std::vector<double>> multipleEulerMaruyama(stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

std::vector<std::vector<double>> multipleEulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

std::vector<double> averageEulerMaruyama(stochasticModel model, int numSimulations);

#endif // __STOCHASTICMETHODS_H__