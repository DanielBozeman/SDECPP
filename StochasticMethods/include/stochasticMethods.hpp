#ifndef __STOCHASTICMETHODS_H__
#define __STOCHASTICMETHODS_H__

#include <vector>
#include <string>

typedef double (*stochastic_function)(double&, double&, std::vector<double>&);

typedef double (*time_function)(double&, double&, std::vector<double>&, std::vector<double>&);

/**
 * @brief Object representing a stochastic differential equation
 * 
 */
class stochasticModel{
    private:

    public:
    stochastic_function alphaFunction;
    stochastic_function betaFunction;

    time_function timeAlpha;
    time_function timeBeta;

    bool isTime;

    double calculateAlpha(double value, double times);
    double calculateBeta(double value, double times);

    double initialValue;
    std::vector<double> timeInterval;
    std::vector<std::vector<double>> parameters;
    std::vector<std::vector<std::vector<double>>> parameterLimits;
    std::vector<std::vector<double>> parameterSteps;

    stochasticModel(stochastic_function function1, stochastic_function function2, double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<std::vector<double>>> constantLimits = {}, std::vector<std::vector<double>> stepSizes = {});  
    
    stochasticModel(time_function function1, time_function function2, double startValue, std::vector<double> times, std::vector<std::vector<double>> constants, std::vector<std::vector<std::vector<double>>> constantLimits = {}, std::vector<std::vector<double>> stepSizes = {});  
    
    void setParameters(std::vector<std::vector<double>> constants);  
    virtual void randomizeParameter(int paramSet);
    virtual void parameterNeighbor(int paramSet);
};  

class polynomialModel : public stochasticModel{
    
    private:

    public:

    std::vector<std::vector<int>> activeTerms;

    polynomialModel(double startValue, std::vector<double> times, std::vector<std::vector<double>> constants = {{},{},{}}, std::vector<std::vector<int>> usedTerms= {{},{},{}}, std::vector<std::vector<std::vector<double>>> constantLimits = {{},{},{}}, std::vector<std::vector<double>> stepSizes = {{},{},{}});  

    void addNextTerm(int paramSet);
    void addTerm(int paramSet, int term);
    void addMultipleTerms(int paramSet, std::vector<int> terms);
    void addMultipleTerms(int paramSet, int maxTerm);


    void removeLastTerm(int paramSet);
    void removeTerm(int paramSet, int term);
    void removeAllTerms(int paramSet);
    void removeAllTerms();


    void setTermParameter(int paramSet, int term, double coefficient);
    void randomizeParameter(int paramSet) override;
    void parameterNeighbor(int paramSet) override;

    std::string toString();

    std::vector<long double> calculateAIC(std::vector<double> &observations, int numSims = 100000, int divisions = 100, double percentage = 0.02);
};

double zeroFunction(double& value, double& time, std::vector<double>& parameters);

double polynomialFunction(double& value, double& time, std::vector<double>& parameters);

double polynomialTimeFunction(double& value, double& time, std::vector<double>& valueParams, std::vector<double>& timeParams);

double zeroTimeFunction(double& value, double& time, std::vector<double>& valueParams, std::vector<double>& timeParams);

double polynomialNoTimeFunction(double&value, double& time, std::vector<double>& parameters, std::vector<double>& timeParams);

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