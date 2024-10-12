#include "RandomUtils.hpp"

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