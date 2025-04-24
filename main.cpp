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

int binarySearch(std::vector<double>& array, double value){
    int low = 0;
    int high = array.size() - 1;
    int rank = 0;

    while(low <= high){
        int mid = low + (high - low)/2;

        if(array[mid] < value){
            rank = mid + 1;
            low = mid + 1;
        }else{
            high = mid - 1;
        }
    }

    return rank;
}

void rewriteTest(){

    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    std::vector<std::vector<double>> parameters = {{1},{0.2}};
    std::vector<std::vector<std::vector<double>>> parameterLimits = {{{0,2}},{{0,1}}};
    std::vector<std::vector<double>> parameterSteps = {{1},{1}};

    double initialValue = 1;

    std::vector<double> times;

    double start = 0;
    double end = 1;
    double dt = 0.005;
    double divisions = 2;

    linearlySpacedVector(times, start, end, dt);

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, parameters, parameterLimits, parameterSteps);

    std::vector<double> output;

    std::vector<double> emOutput;

    eulerMaruyamaByReference(output, model);

    model.parameters[0] = {0};

    model.parameters[0] = paramEstimation(model, 0, output, 1, 200, 0.99, 250, 15, driftCost);

    std::cout << "\nDrift: " << model.parameters[0][0];

    model.parameters[1] = {0.0001};

    model.parameters[1] = paramEstimation(model, 1, output, 1, 100, 0.9, 100, 1, varianceCost, {100, 10});

    std::cout << "\nParam: " << model.parameters[1][0];


}

void testLikelihood(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto fakeAlpha = [](double& value, double& time, std::vector<double>& parameters){
        return(parameters[0] * value + parameters[1]);
    };

    std::vector<std::vector<double>> trueParameters = {{1},{0.2}};
    std::vector<std::vector<std::vector<double>>> parameterLimits = {{{0,2}},{{0,1}}};
    std::vector<std::vector<double>> parameterSteps = {{1},{1}};



    std::vector<std::vector<double>> falseParameters = {{0,0},{0.001}};
    std::vector<std::vector<std::vector<double>>> falseParameterLimits = {{{0,2},{-100,100}},{{0,1}}};
    std::vector<std::vector<double>> falseParameterSteps = {{1,1},{1}};


    double initialValue = 1;

    std::vector<double> times;

    double start = 0;
    double end = 5;
    double dt = 0.05;
    double divisions = 2;

    linearlySpacedVector(times, start, end, dt);

    stochasticModel trueModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, trueParameters, parameterLimits, parameterSteps);

    stochasticModel falseModel = stochasticModel(fakeAlpha, betaFunction, initialValue, times, falseParameters, falseParameterLimits, falseParameterSteps);

    std::vector<double> output;

    eulerMaruyamaByReference(output, trueModel);

    for(int i = 0; i < output.size(); i++){
        std::cout << "\n" << output[i];
    }

    trueModel.parameters = {{0},{0.0001}};

    trueModel.parameters[0] = paramEstimation(trueModel, 0, output, 1, 200, 0.99, 250, 15, driftCost);

    std::cout << "\n===============================";

    trueModel.parameters[1] = paramEstimation(trueModel, 1, output, 1, 100, 0.9, 100, 1, varianceCost, {20, 100});

    std::cout << "\n===============================";


    falseModel.parameters[0] = paramEstimation(falseModel, 0, output, 1, 200, 0.99, 250, 15, driftCost);

    std::cout << "\n===============================";


    falseModel.parameters[1] = paramEstimation(falseModel, 1, output, 1, 100, 0.9, 100, 1, varianceCost, {20, 100});

    std::cout << "\n===============================";


    std::cout << "\nFinding Likelihood";

    std::vector<std::vector<double>> newParams = {{0.5},{0.2}};
    stochasticModel newModel = trueModel;
    newModel.parameters = newParams;


}

void pdfEstimation(){
    std::vector<double> standardNormalVariables = {};

    for(int i = 0; i < 100000; i++){
        standardNormalVariables.push_back(randomPathMaker::dW(1));
    }

    double input = 4;

    double estimate = estimatePdf(standardNormalVariables, input, 0.02);
    //double estimate = 1;
    double truePdf = normalPDF(input, 0, 1);
    
    std::cout << "\nEstimate: " << estimate;
    std::cout << "\nPDF: " << truePdf;
}

void testPath(){
    std::vector<double> standardNormals = {};

    for(int i = 0; i < 100000; i++){
        standardNormals.push_back(randomPathMaker::dW(1));
    }

    double sum = 0;
    for(int i = 0; i < standardNormals.size(); i++){
        sum += standardNormals[i];
    }

    double mean = sum/standardNormals.size();

    std::cout << "\nMean is: " << mean;

    std::sort(standardNormals.begin(), standardNormals.end());

    std::cout << "\nMin: " << standardNormals[0] << "\nMax: " << standardNormals.back();
}

void testFits(){
   double input = -9.8;

   long double mantissa = 0;
   int exp = 0;

   mantissa = frexp(-9.8, &exp);

   long unsigned int exponent = 0;

   exponent = exp;

   std::cout << "\nMantissa: " << mantissa << "   Exp: " << exponent;


}

void testPolynomial(){

    double initialValue = 1;

    double start = 0;
    double end = 1;
    double dt = 0.05;
    int divisions = 10;

    std::vector<double> times;

    linearlySpacedVector(times, start, end, dt);
    
    polynomialModel polyModel = polynomialModel(polynomialWithXFunction, polynomialNoTimeFunction,initialValue, times);

    polyModel.addTerm(0, 0);

    polyModel.setTermParameter(0, 0, 1);


    std::vector<double> testOutput;
    std::vector<double> output;

    std::cout << "\n" << polyModel.toString();

    eulerMaruyamaWithin(output, polyModel, divisions);

    for(int i = 0; i < output.size(); i++){
        std::cout << "\n" << output[i];
    }

    polynomialModel testModel = polynomialModel(polynomialWithXFunction, polynomialNoTimeFunction, output[0], times);

    std::cout << "\nFitting!";

    testModel = bestModelNTerms(output, times, 1, 5, divisions);

    eulerMaruyamaWithin(testOutput, testModel, divisions);

    // std::vector<long double> testAIC = polyModel.calculateAIC(output, 25000, divisions);
    // std::vector<long double> aic2 = testModel.calculateAIC(output, 25000, divisions);

    // std::cout << "\nReal AIC: " << testAIC[0];
    // std::cout << "\nFit AIC: " << aic2[0];

    std::cout << "\n" << testModel.toString();
    //polynomialModel model = fitPolynomial(output, times, 5, divisions);
    
    std::vector<std::vector<double>> approximations;

    multipleEulerMaruyamaWithin(approximations, testModel, divisions, 1000);

    multiVectorToCSV(approximations, "estimationOutput.csv");

    vectorToCSV(output, "dataOutput.csv");

    std::cout << "\nFinished";

}

void constructionTests(){
    polynomialModel testModel = polynomialModel(polynomialWithXFunction, polynomialNoTimeFunction, 0, {0,1,2});

    testModel.addMultipleRandomTerms({0,2}, 5, 5);
    //testModel.addRandomTerm(2, 5);

    std::cout << "\n" << testModel.toString();

    // for(int i = 0; i < 6; i++){
    //     testModel.neighboringSet(5);

    //     std::cout << "\nNew Model: \n" << testModel.toString();
    // }
}

void bugTesting(){

    double initialValue = 1;

    double start = 0;
    double end = 1;
    double dt = 0.05;

    std::vector<double> times;

    linearlySpacedVector(times, start, end, dt);

    polynomialModel polyModel = polynomialModel(initialValue, times);

    polyModel.addTerm(0, 1);
    polyModel.addTerm(1, 1);

    polyModel.setTermParameter(0, 1, 1);
    polyModel.setTermParameter(1, 1, 0.17);

    std::vector<double> output;

    eulerMaruyamaWithin(output, polyModel, 100);

    std::cout << "\nFinished!";   
}

void infinityTest(){
    double inf = INFINITY;
    inf = std::numeric_limits<double>::max();
    double time = 0.05;

    std::vector<double> params = {0,0,0,0.05};

    double output = polynomialFunction(inf, time, params);

    std::cout << "\nDone";
}

int main(){

    #ifdef __OPTIMIZE__
    std::cout << "\nOptimization enabled!";
    #else
    std::cout << "\nNo optimization set!";
    #endif
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

    //rewriteTest();

    //testLikelihood();

    //testPdf();
    //pdfEstimation();

    //testFits();

    testPolynomial();

    //constructionTests();

    //bugTesting();

    //infinityTest();

    //testPath();

    std::cout << "\nDone";
}