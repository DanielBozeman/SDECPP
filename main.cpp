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
    double variance = 100;

    std::vector<double> dWs = {};

    for(int i = 0; i < 100000; i++){
        dWs.push_back(randomPathMaker::dW(variance));
    }

    std::cout << "\nVar is " << sampleVariance(dWs) << "    mean is " << sampleMean(dWs);
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
    double intervalEnd = 1;

    double timeDiscretization = 0.003;

    std::vector<std::vector<double>> constants = {{0.15}, {0.3}};

    double initialValue = 120;

    std::vector<double> times = {0};

    times.reserve((intervalEnd - intervalStart)/timeDiscretization);

    while (times.back() < intervalEnd){
        times.push_back(times.back()+timeDiscretization);
    }

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants);

    std::vector<std::vector<double>> approximation = multipleEulerMaruyama(model, numSims);

    multiVectorToCSV(approximation, "output.csv");

    return {};
}

void fileWriterTest(){

    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<std::vector<double>> stocks = {stockCloses};

    multiVectorToCSV(stocks, "output.csv");
}

void rmseTest(){
    std::vector<std::vector<double>> actual = {{0,4,5,6},{0,2,5,6}};
    std::vector<double> prediction = {0,0,0,0};

    std::cout << multiVectorRMSE(actual, prediction);

}

std::vector<std::vector<double>> simulatedAnnealingTest(){
    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);

    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());

    stockCloses = y;

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    // for(int i = 0; i < stockCloses.size(); i++){
    //     std::cout << "\nStock Value: " << stockCloses[i];
    // }

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

    std::vector<double> trueReturns = {};

    for(int i = 1; i < observations.size(); i++){
        trueReturns.push_back((observations[i] - observations[i-1]));
    }

    while(simModel.parameters[1][0] < 1){
        std::vector<std::vector<double>> sims = multipleEulerMaruyama(simModel, 500);

        double cost = returnComparison(sims, trueReturns);

        std::cout << "\nVariance: " << simModel.parameters[1][0] <<  "    Cost: " << cost;

        simModel.parameters[1][0] += 0.0002;
    }
}

void stats(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * pow(value, parameters[1]));
    };

    std::vector<double> times = {};
    times.push_back(0);

    while(times.back() < 20){
        times.push_back(times.back() + 0.01);
    }

    double initialValue = 1.0;

    std::vector<std::vector<double>> baseConstants = {{0.05},{0.15,0.5}};

    stochasticModel obsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    std::vector<double> observations = eulerMaruyama(obsModel);

    stochasticModel noVarModel = obsModel;
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

    std::vector<std::vector<double>> constants = {{0.05}, {0.08, 0.5}};

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants);

    for(int i = 0; i < 100; i++){
        std::vector<std::vector<double>> curApproximation = multipleEulerMaruyama(model, 500);

        for(int j = 0; j < curApproximation.size(); j++){
            for(int k = 0; k < curApproximation[k].size(); k++){
                curApproximation[j][k] = curApproximation[j][k] - drift[k];
            }
        }

        double RMS = returnComparison(curApproximation, trueReturns);

        std::cout << "\nCurrent Var: " << "  " << model.parameters[1][0] << "    " << model.parameters[1][1] << "    Cost: " << RMS;
        
        model.parameters[1][0] += 0.001;
    }
    return;
}

void cevTest(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * pow(value, parameters[1]));
    };

    std::vector<double> times = {};
    times.push_back(0);

    while(times.back() < 20){
        times.push_back(times.back() + 0.01);
    }

    double initialValue = 1.0;

    std::vector<std::vector<double>> baseConstants = {{0.05},{0.15,0.5}};

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    stochasticModel noVarModel = model;
    noVarModel.betaFunction = zeroFunction;

    std::vector<double> drift = eulerMaruyama(noVarModel);

    std::vector<double> observations = eulerMaruyama(model);

    for(int i = 0; i < observations.size(); i++){
        observations[i] = observations[i] - drift[i];
    }

    std::vector<double> trueReturns = {};

    for(int i = 1; i < observations.size(); i++){
        trueReturns.push_back((observations[i] - observations[i-1]));
    }


    std::vector<std::vector<double>> simulations = multipleEulerMaruyama(model, 1000);

    for(int i = 0; i < simulations.size(); i++){
        for(int j = 0; j < simulations[i].size(); j++){
            simulations[i][j] = simulations[i][j] - drift[j];
        }
    }

    std::vector<std::vector<double>> simReturns = {};

    for(int i = 1; i < simulations[0].size(); i++){
        simReturns.push_back({});
        for(int j = 0; j < simulations.size(); j++){
            simReturns[i-1].push_back(simulations[j][i] - simulations[j][i-1]);
        }
    }

    

    //std::cout << "\nExample sim 9-10 " << simulations[1][9] << " to " << simulations[1][10] << " and real return is " << observations[9] << " to " << observations[10];

    //std::cout << "\nMean is " << sampleMean(simReturns[1999]);

    for(int i = 0; i < trueReturns.size(); i++){
        double simMean = sampleMean(simReturns[i]);
        double simVariance = sampleVariance(simReturns[i]);
        std::cout << "\ni = " << i;
        double pdf = normalPDF(trueReturns[i], simMean, simVariance);
        std::cout << "\nPDF is " << pdf;
        std::cout << "\nNormal is mean: " << simMean << " and var: " << simVariance << "|| Observed return is: " << trueReturns[i];
    }
}

void testing(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * pow(value, parameters[1]));
    };

    std::vector<double> times = {};
    times.push_back(0);

    while(times.back() < 12){
        times.push_back(times.back() + 0.01);
    }

    double initialValue = 30.0;

    std::vector<std::vector<double>> baseConstants = {{0.05},{1.5, 0.4}};

    stochasticModel obsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    std::vector<double> stockCloses = eulerMaruyama(obsModel);

    std::vector<std::vector<double>> approx = multipleEulerMaruyama(obsModel, 500);

    vectorToCSV(stockCloses, "dataOutput.csv");

    std::vector<double> trueReturns = {};

    for(int i = 1; i < stockCloses.size(); i++){
        trueReturns.push_back((stockCloses[i] - stockCloses[i-1]));
    }

    multiVectorToCSV(approx, "estimationOutput.csv");

    std::cout << "\nCost: " << returnComparison(approx, trueReturns);
    std::cout << "\nSize: " << stockCloses.size();
    std::cout << "\nTimes: " << times.size();
}

void testBS(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    std::vector<std::vector<double>> constants = {{0.000891396},{0.00951891}};

    std::vector<double> stockCloses = csvColumnToVector("StockData/SPX_Post61.csv", 6);
    
    std::vector<double> y(stockCloses.end() - 500, stockCloses.end());
    stockCloses = y;

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    double initialValue = stockCloses[0];

    stochasticModel bsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants);

    bsModel.parameters[1][0] /= 2;
    
    double stat = testFit(bsModel, stockCloses, 90);

    std::cout << "\nParameter is " << bsModel.parameters[1][0] << "   stat is " << stat;

    bsModel.parameters[1][0] *= 2;

    stat = testFit(bsModel, stockCloses, 90);

    std::cout << "\nParameter is " << bsModel.parameters[1][0] << "   stat is " << stat;

    bsModel.parameters[1][0] *= 2;

    stat = testFit(bsModel, stockCloses, 90);

    std::cout << "\nParameter is " << bsModel.parameters[1][0] << "   stat is " << stat;
}

//--------------------------------------------

void rewriteTest(){

    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    std::vector<std::vector<double>> parameters = {{1},{0.2}};

    double initialValue = 1;

    std::vector<double> times = {};

    double start = 0;
    double end = 10;
    double dt = 1;
    
    linearlySpacedVector(times, start, end, dt);

    std::vector<double> output;

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, parameters);
    
    randomPathMaker rp = randomPathMaker();

    std::vector<double> path = rp.makePath(times[0], times.back() + dt, dt);

    std::vector<std::vector<double>> pathSet = {};

    double sum = 0;

    int numSims = 1000000;

    path = rp.makePath(times[0], times.back() + dt, dt); 

    for(int i = 0; i < numSims; i++){
        path = rp.makePath(times[0], times.back()+ dt, dt);
        sum += path.back();
    }

    sum = sum / numSims;

    std::cout << "\nAverage dW: " << sum;

    //multiVectorToCSV(pathSet, "output.csv");
    return;
}

int main(){
    //SAComparison();
    //varianceViewer();
    //logDifferenceTest();
    //statsTests();
    //driftVolFinder();
    //stats();
    //fitOrnstein("StockData/SPX_Post61.csv", 6, -500);
    //fitBlackScholes("StockData/SPX_Post61.csv", 6, -500);
    //rmseTest();^
    //eulerMaruyamaTest(500);
    //createPath();
    //brownianPathTester();
    //fileWriterTest();
    //testing();
    //fitCEV();
    //cevTest();
    //fitRandom("StockData/SPX_Post61.csv", 6, -500);
    //testBS();

    rewriteTest();
}