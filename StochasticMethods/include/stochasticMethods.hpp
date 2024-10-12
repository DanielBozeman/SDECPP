#ifndef __STOCHASTICMETHODS_H__
#define __STOCHASTICMETHODS_H__

#include <vector>

typedef double (*stochastic_function)(double, double, std::vector<double>);

class stochasticModel{
    private:

    public:
    stochasticModel(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters);
    stochastic_function alphaFunction;
    stochastic_function betaFunction;
    double initialValue;
    std::vector<double> timeInterval;
    std::vector<double> parameters;
};

std::vector<double> eulerMaruyama(stochasticModel model, std::vector<double> brownianPath);

std::vector<double> eulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, std::vector<double> brownianPath = {});

std::vector<std::vector<double>> multipleEulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<std::vector<double>> brownianPaths, std::vector<double> parameters, int numSimulations);

#endif // __STOCHASTICMETHODS_H__