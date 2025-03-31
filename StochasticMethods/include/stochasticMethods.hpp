#ifndef __STOCHASTICMETHODS_H__
#define __STOCHASTICMETHODS_H__

#include <vector>

typedef double (*stochastic_function)(double&, double&, std::vector<double>&);

/**
 * @brief Object representing a stochastic differential equation
 * 
 */
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

    /**
     * @brief Construct a new stochastic Model object
     * 
     * @param function1 Alpha function for the SDE
     * @param function2 Beta function of the SDE
     * @param startValue Initial value of the SDE
     * @param times Vector representing each timestep that the SDE will be simulated over
     * @param constants Vector of vectors, first vector is the parameters for the alpha function and second is the beta function parameters
     * @param constantLimits Upper and lower bounds for each parameter
     * @param stepSizes How far apart neighbors are in the simulated annealing
     */
    stochasticModel(stochastic_function function1, stochastic_function function2, double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<std::vector<double>>> constantLimits = {}, std::vector<std::vector<double>> stepSizes = {});  
    void setParameters(std::vector<std::vector<double>> constants);  
};

double zeroFunction(double& value, double& time, std::vector<double>& parameters);

void linearlySpacedVector(std::vector<double> &xs, double a, double b, double h);

void linearlySpacedVectorBySize(std::vector<double> &xs, double a, double b, std::size_t N);

std::vector<double> eulerMaruyama(stochasticModel model, std::vector<double> brownianPath = {});

std::vector<double> eulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, std::vector<double> brownianPath = {});

void eulerMaruyamaWithin(std::vector<double> &approximation, stochasticModel model, int numDivisions, std::vector<double> brownianPath = {});

void eulerMaruyamaByReference(std::vector<double> &approximation, stochasticModel model, std::vector<double> brownianPath = {});

void multipleEulerMaruyamaWithin(std::vector<std::vector<double>> &approximations, stochasticModel model, int numDivisions, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

std::vector<std::vector<double>> multipleEulerMaruyama(stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

std::vector<std::vector<double>> multipleEulerMaruyama(stochastic_function alphaFunction, stochastic_function betaFunction, double initialValue, std::vector<double> timeInterval, std::vector<double> parameters, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

void multipleEulerMaruyamaByReference(std::vector<std::vector<double>> &approximations, stochasticModel model, int numSimulations, std::vector<std::vector<double>> brownianPaths = {});

std::vector<double> averageEulerMaruyama(stochasticModel model, int numSimulations);

double testFit(stochasticModel model, std::vector<double> data, int numInterRuns);

#endif // __STOCHASTICMETHODS_H__