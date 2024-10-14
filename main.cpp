#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
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

double alphaFunction(double value, double time, std::vector<double> parameters){
    return (parameters[0]*value);
}

double betaFunction(double value, double time, std::vector<double> parameters){
    return (parameters[1]*value);
}

std::vector<double> eulerMaruyamaTest(int numSims){

    randomPathMaker rp;

    double intervalStart = 0;
    double intervalEnd = 10;

    double timeDiscretization = pow(2, -8);

    std::vector<double> constants = {0.05, 0.2};

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

    double rms = rmse(actual, prediction);

    std::cout << "\n RMSE is " << rms;
}

void parameterNeighborTest(){
    std::vector<double> parameters = {0, 5, 7};
    std::vector<std::vector<double>> parameterLimits = {{-2,0},{0,20},{-20,8}};
    std::vector<double> parameterSteps = {0.1, 1, 0.001};

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "\nStarting test";

    for(int i = 0; i < 1000; i++){
        //std::cout << "\n";
        //for(int j = 0; j < parameters.size(); j++){
            //std::cout << parameters[j] << " ";
        //}
        parameterNeighbor(parameters, parameterLimits, parameterSteps);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "\nTime taken: " << ms_double.count()/1000;
}

void simulatedAnnealingTest(){
    std::vector<double> observations = {0,2,6,8,3,15,8,3,7,0,2};

    randomPathMaker rp;

    double intervalStart = 0;
    double intervalEnd = 10;

    double timeDiscretization = 1;

    std::vector<double> constants = {0.05, 0.2};

    std::vector<std::vector<double>> constantLimits = {{-4,5},{0,1}};

    std::vector<double> constantSteps = {0.1, 0.001};

    double initialValue = 5;

    std::vector<double> times = {0};

    times.reserve((intervalEnd - intervalStart)/timeDiscretization);

    while (times.back() < intervalEnd){
        times.push_back(times.back()+timeDiscretization);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    simulatedAnnealingParameterEstimation(model, observations, 500, 100, 0.9, 50, 0.1);

    std::cout << "\n\nParameter 1: " << model.parameters[0];
    std::cout << "\nParameter 2: " << model.parameters[1]; 
}

int main(){
    std::vector<double> stockRows = csvColumnToVector("StockData/SPX_Post61.csv", 0);

    for(int i = 0; i < stockRows.size(); i ++){
        std::cout << "\n" << stockRows[i];
    }
}