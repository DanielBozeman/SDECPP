#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include "RandomUtils.hpp"
#include "stochastic.hpp"
#include "FileUtils.hpp"

void brownianPathTester(){
    randomPathMaker rp;

    std::string fileName = "output.csv";

    int status = remove("output.csv");

    std::ofstream outfile(fileName);

    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> paths = rp.makeMultiplePaths(0, 10, 0.001, 500);

    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "\nTime: " << ms_double.count()/1000;
}

double alphaFunction(double& value, double& time, std::vector<double>& parameters){
    return (parameters[0]*value);
}

double betaFunction(double& value, double& time, std::vector<double>& parameters){
    return (parameters[0]*value);
}

std::vector<double> eulerMaruyamaTest(int numSims){

    randomPathMaker rp;

    double intervalStart = 0;
    double intervalEnd = 10;

    double timeDiscretization = pow(2, -8);

    std::vector<std::vector<double>> constants = {{0.05}, {0.2}};

    double initialValue = 5;

    std::vector<double> times = {0};

    times.reserve((intervalEnd - intervalStart)/timeDiscretization);

    while (times.back() < intervalEnd){
        times.push_back(times.back()+timeDiscretization);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants);

    std::vector<double> approximation = averageEulerMaruyama(model, numSims);

    return approximation;
}

void fileWriterTest(){

    std::vector<std::vector<double>> vectors = {};

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "\nStarting test";

    std::vector<double> approximation = eulerMaruyamaTest(1000);

    vectorToCSV(approximation, "output.csv");

    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "\nTime taken: " << ms_double.count()/1000;
}

void rmseTest(){
    std::vector<double> actual = {0,3,5,6};
    std::vector<double> prediction = {0,0,0,0};

}

std::vector<std::vector<double>> simulatedAnnealingTest(){
    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    std::vector<std::vector<double>> constants = {{0}, {0.0001}};

    std::vector<std::vector<std::vector<double>>> constantLimits = {{{0,2}},{{0.0001,2}}};

    std::vector<std::vector<double>> constantSteps = {{0.01}, {0.01}};

    double initialValue = 3655.04;

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    for(int i = 0; i < 3; i++){

        std::cout << "\nStarting Drift Estimation";
        model.parameters[0] = simulatedAnnealingDriftEstimation(model, 0, stockCloses, 1, 200, 0.9, 150, 1, multiVectorRMSE);
        
        //model.parameters[0] = {0.000865574};
        std::cout << "\nDrift Est: " << model.parameters[0][0];
        std::cout << "\nStarting Vol Estimation";

        model.betaFunction = betaFunction;
        model.parameters[1] = simulatedAnnealingVolEstimation(model, 1, stockCloses, 500, 20, 0.9, 150, 10, returnComparison);
        
        std::cout << "\n\n For Run " << i;
        std::cout << "\nDrift Est: " << model.parameters[0][0];
        std::cout << "\nVolatility Est: " << model.parameters[1][0]; 

        for(int j = 0; j < constantSteps.size(); j++){
            for(int k = 0; k < constantSteps[j].size(); k++){
                constantSteps[j][k] *= 0.01;
            }
        }

        model.parameterSteps = constantSteps;
    }
    

    return model.parameters;
}

std::vector<std::vector<double>> driftVolFinder(){
    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<double> returns;

    for (int i = 1; i < stockCloses.size(); i++){
        double baseReturn = stockCloses[i]/stockCloses[i-1];
        returns.push_back(log(baseReturn));
    }

    double average = std::accumulate(returns.begin(), returns.end(),0.0) / returns.size();

    double variance = 0;

    for(int i = 0; i < returns.size(); i++){
        variance += ((returns[i] - average)*(returns[i]-average));
    }
    variance /= returns.size();

    variance = sqrt(variance);    

    double logAverage = std::accumulate(returns.begin(), returns.end(),0.0) / returns.size();

    std::cout << "\nVolatility: " << variance;
    std::cout << "\nDrift: " << logAverage;

    std::vector<std::vector<double>> parameters = {{logAverage}, {variance}};
    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, 3655.04, times, parameters);

    std::vector<double> approximation = averageEulerMaruyama(model, 1500);

    std::vector<std::vector<double>> lines = {approximation, stockCloses};

    return lines;
}

void varianceViewer(){
    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<double> returns;

    for (int i = 1; i < stockCloses.size(); i++){
        double baseReturn = stockCloses[i]/stockCloses[i-1];
        returns.push_back(log(baseReturn));
    }

    double average = std::accumulate(returns.begin(), returns.end(),0.0) / returns.size();

    double variance = 0;

    for(int i = 0; i < returns.size(); i++){
        variance += ((returns[i] - average)*(returns[i]-average));
    }
    variance /= returns.size();

    variance = sqrt(variance);   

    double logAverage = std::accumulate(returns.begin(), returns.end(),0.0) / returns.size();

    std::cout << "\nVolatility: " << variance;
    std::cout << "\nDrift: " << logAverage;

    std::vector<std::vector<double>> parameters = {{logAverage}, {variance}};
    //std::vector<std::vector<double>> parameters = {{logAverage}, {0.005}};
    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, 3655.04, times, parameters);

    std::vector<std::vector<double>> trajectories = multipleEulerMaruyama(model, 1000);


    trajectories.push_back(stockCloses);

    multiVectorToCSV(trajectories, "output.csv");

    return;
}

void SAComparison(){
    std::vector<std::vector<double>> parameters = simulatedAnnealingTest();
    std::vector<std::vector<double>> lines = driftVolFinder();

    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, 3655.04, times, parameters);

    std::vector<double> newApproximation = averageEulerMaruyama(model, 500);

    lines.push_back(newApproximation);

    multiVectorToCSV(lines, "output.csv");
}

void statsTests(){

    double simVariance = 0.005;
    double trueVariance = 0.00951891;
    double startValue = 3655.04;
    double drift = 0.000891396;

    int timeLength = 500;
    std::vector<double> times = {};

    for(int i = 0; i < timeLength; i++){
        times.push_back(i);
    }

    std::vector<std::vector<double>> simConstants = {{drift},{simVariance}};
    std::vector<std::vector<double>> observationConstants = {{drift},{trueVariance}};

    stochasticModel simModel = stochasticModel(zeroFunction, betaFunction, startValue, times, simConstants);
    stochasticModel obsModel = stochasticModel(zeroFunction, betaFunction, startValue, times, observationConstants);
    
    std::vector<double> observations = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(observations.end() - 500, observations.end());

    observations = y;

    std::vector<std::vector<double>> costs = {};

    std::vector<std::vector<double>> outReturns = {};

    for(int i = 0; i < 50; i++){
        std::vector<std::vector<double>> sims = multipleEulerMaruyama(simModel, 500);

        double cost = returnComparison(sims, observations);

        std::cout << "\nVariance: " << simModel.parameters[1][0] <<  "    Cost: " << cost;

        simModel.parameters[1][0] += 0.0002;
    }

    //multiVectorToCSV(outReturns, "output.csv");

    //multiVectorToCSV(sims, "output.csv");

    //vectorToCSV(sims[10], "output.csv");


}

void stats(){
    std::vector<double> pdfTests = {2.5, 5, 10, 25};

    double mean = 12;
    double variance = 25;

    for(int i = 0; i < pdfTests.size(); i++){
        std::cout << "\nPDF of " << pdfTests[i] << " is " << normalPDF(pdfTests[i], mean, variance);
    }
}

int main(){
    SAComparison();
    //varianceViewer();
    //logDifferenceTest();
    //statsTests();
    //driftVolFinder();
    //stats();
}