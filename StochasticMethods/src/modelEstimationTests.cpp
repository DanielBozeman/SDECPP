#include "modelEstimationTests.hpp"
#include "parameterEstimation.hpp"
#include "stochasticMethods.hpp"
#include "FileUtils.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>


stochasticModel fitBlackScholes(std::string fileName, int dataColumn, int dataStart, int dataEnd){

    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0]*value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0]*value);
    };
    
    std::vector<double> stockCloses = csvColumnToVector(fileName, dataColumn);

    std::vector<double> y = {};
    
    if(dataStart < 0){
        std::vector<double> y(stockCloses.end() + dataStart, stockCloses.end() - dataEnd);
        stockCloses = y;
    }else if(dataEnd < 0){
        std::vector<double> y(stockCloses.begin() + dataStart, stockCloses.begin() - dataEnd);
        stockCloses = y;
    }else{
        std::vector<double> y(stockCloses.begin() + dataStart, stockCloses.end() - dataEnd);
        stockCloses = y;
    }

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    double initialValue = stockCloses[0];

    std::vector<std::vector<double>> constants = {{0}, {0.0001}};

    //Starting with a very wide constant limit and step values
    std::vector<std::vector<std::vector<double>>> constantLimits = {{{0,10}},{{0.0001,10}}};

    std::vector<std::vector<double>> constantSteps = {{1}, {1}};

    std::cout << "\nInitial Value: " << initialValue;

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    stochasticModel noVarModel = stochasticModel(alphaFunction, zeroFunction, initialValue, times, constants, constantLimits, constantSteps);

    for(int i = 0; i < 5; i++){

        std::cout << "\nStarting Drift Estimation";
      
        //model.parameters[0] = paramEstimation(model, 0, stockCloses, 1, 2000, 0.9, 150, 100, driftCost);
        
        model.parameters[0] = {0.000865574};
        //std::cout << "\nDrift Est: " << model.parameters[0][0];
        std::cout << "\nStarting Vol Estimation";

        model.betaFunction = betaFunction;
        //model.parameters[1] = simulatedAnnealingVolEstimation(model, 1, stockCloses, 500, 20, 0.9, 200, 10, returnComparison);
        model.parameters[1] = paramEstimation(model, 1, stockCloses, 500, 20, 0.9, 200, 10, varianceCost);

        //std::cout << "\n\n For Run " << i;
        //std::cout << "\nDrift Est: " << model.parameters[0][0];
        //std::cout << "\nVolatility Est: " << model.parameters[1][0]; 

        for(int j = 0; j < constantSteps.size(); j++){
            for(int k = 0; k < constantSteps[j].size(); k++){
                constantSteps[j][k] *= 0.01;
            }
        }

        model.parameterSteps = constantSteps;
    }

    //double stat = testFit(model, stockCloses, 90);

    //std::cout << "\nTest stat is " << stat;

    multiVectorToCSV({stockCloses}, "dataOutput.csv");

    std::vector<std::vector<double>> estimationVectors = multipleEulerMaruyama(model, 500);

    multiVectorToCSV(estimationVectors, "estimationOutput.csv");

    std::vector<std::vector<double>> trueConstants = {{0.000891396}, {0.00951891}};

    stochasticModel trueModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, trueConstants);

    std::vector<std::vector<double>> trueVectors = multipleEulerMaruyama(trueModel, 500);

    multiVectorToCSV(trueVectors, "trueOutput.csv");

    return model;
}

stochasticModel fitOrnstein(std::string fileName, int dataColumn, int dataStart, int dataEnd){
    //Ornstein is dx_t = theta * (mu - x_t) dt +  sigma * dW_t

    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0]*(parameters[1] - value));
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0]);
    };

    //Constants are ordered as {theta, mu},{sigmat}
    std::vector<std::vector<double>> constants = {{0.0001,0.0001}, {0.0001}};

    //Starting with a very wide constant limit and step values
    std::vector<std::vector<std::vector<double>>> constantLimits = {{{0.0001,10},{-10,10}},{{0.0001,10}}};

    std::vector<std::vector<double>> constantSteps = {{1,1}, {1}};

    std::vector<double> times = {};
    times.push_back(0);

    // for(int i = 0; i < 12; i++){
    //     times.push_back(i);
    // }

    while(times.back() < 12){
        times.push_back(times.back() + 0.01);
    }

    //double initialValue = stockCloses[0];
    double initialValue = 10.0;

    std::vector<std::vector<double>> baseConstants = {{0.8,1.0},{1.44}};

    stochasticModel obsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    std::vector<double> stockCloses = eulerMaruyama(obsModel);

    vectorToCSV(stockCloses, "dataOutput.csv");

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    stochasticModel baseModel = model;

    for(int i = 0; i < 1; i++){

        std::cout << "\nStarting Drift Estimation";
        //model.parameters[0] = simulatedAnnealingDriftEstimation(model, 0, stockCloses, 1, 200, 0.9, 150, 1, multiVectorRMSE);
        model.parameters[0] = paramEstimation(model, 0, stockCloses, 1, 200, 0.9, 150, 1, driftCost);

        //model.parameters[0] = {0.000865574};
        std::cout << "\nDrift Est: " << model.parameters[0][0] << "   " << model.parameters[0][1];

        std::cout << "\nStarting Vol Estimation";

        //model.betaFunction = betaFunction;
        //model.parameters[1] = simulatedAnnealingVolEstimation(model, 1, stockCloses, 500, 20, 0.9, 150, 14, returnComparison);
        model.parameters[1] = paramEstimation(model, 1, stockCloses, 500, 20, 0.9, 150, 14, varianceCost, {100, 10});
        
        std::cout << "\n\n For Run " << i;
        std::cout << "\nDrift Est: " << model.parameters[0][0] << " " << model.parameters[0][1];
        std::cout << "\nVolatility Est: " << model.parameters[1][0]; 

        for(int j = 0; j < constantSteps.size(); j++){
            for(int k = 0; k < constantSteps[j].size(); k++){
                constantSteps[j][k] *= 0.01;
            }
        }

        model.parameterSteps = constantSteps;
    }

    std::vector<std::vector<double>> estimations = multipleEulerMaruyama(model, 500);

    multiVectorToCSV(estimations, "estimationOutput.csv");

    std::vector<std::vector<double>> trueVectors = multipleEulerMaruyama(obsModel, 500);

    multiVectorToCSV(trueVectors, "trueOutput.csv");

    return model;
}

stochasticModel fitCEV(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * pow(value, parameters[1]));
    };

    //Constants are ordered as {theta, mu},{sigmat}
    std::vector<std::vector<double>> constants = {{0.0001}, {0.0001,0.0001}};

    //Starting with a very wide constant limit and step values
    std::vector<std::vector<std::vector<double>>> constantLimits = {{{0.0001,10}},{{0.0001,10},{-4,4}}};

    std::vector<std::vector<double>> constantSteps = {{1}, {1,1}};

    std::vector<double> times = {};
    times.push_back(0);

    while(times.back() < 20){
        times.push_back(times.back() + 0.01);
    }

    //double initialValue = stockCloses[0];
    double initialValue = 1.0;

    std::vector<std::vector<double>> baseConstants = {{0.05},{0.15,0.5}};

    stochasticModel obsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    std::vector<double> stockCloses = eulerMaruyama(obsModel);

    std::vector<double> stockClosesToWrite = stockCloses;

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    stochasticModel baseModel = model;

    for(int i = 0; i < 3; i++){

        std::cout << "\nStarting Drift Estimation";
        model.parameters[0] = simulatedAnnealingDriftEstimation(model, 0, stockCloses, 1, 200, 0.9, 400, 1, multiVectorRMSE);
        
        //model.parameters[0] = {0.000865574};
        std::cout << "\nDrift Est: " << model.parameters[0][0];
        std::cout << "\nStarting Vol Estimation";

        model.betaFunction = betaFunction;
        model.parameters[1] = simulatedAnnealingVolEstimation(model, 1, stockCloses, 500, 50, 0.9, 100, 15, returnComparison);
        
        std::cout << "\n\n For Run " << i;
        std::cout << "\nDrift Est: " << model.parameters[0][0] ;
        std::cout << "\nVolatility Est: " << model.parameters[1][0]<< " " << model.parameters[1][1]; 

        for(int j = 0; j < constantSteps.size(); j++){
            for(int k = 0; k < constantSteps[j].size(); k++){
                constantSteps[j][k] *= 0.01;
            }
        }

        model.parameterSteps = constantSteps;
    }

    std::vector<std::vector<double>> estimations = multipleEulerMaruyama(model, 500);

    multiVectorToCSV(estimations, "estimationOutput.csv");

    std::vector<std::vector<double>> trueVectors = multipleEulerMaruyama(obsModel, 500);

    multiVectorToCSV(trueVectors, "trueOutput.csv");

    vectorToCSV(stockClosesToWrite, "dataOutput.csv");

    return model;
}

stochasticModel fitRandom(std::string fileName, int dataColumn, int dataStart, int dataEnd){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value + pow(value, parameters[1]) + parameters[2]);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value + pow(value, parameters[1]) + parameters[2]);
    };

    std::vector<double> stockCloses = csvColumnToVector(fileName, dataColumn);

    std::vector<double> y = {};
    
    if(dataStart < 0){
        std::vector<double> y(stockCloses.end() + dataStart, stockCloses.end() - dataEnd);
        stockCloses = y;
    }else if(dataEnd < 0){
        std::vector<double> y(stockCloses.begin() + dataStart, stockCloses.begin() - dataEnd);
        stockCloses = y;
    }else{
        std::vector<double> y(stockCloses.begin() + dataStart, stockCloses.end() - dataEnd);
        stockCloses = y;
    }

    std::vector<double> times = {};

    for(int i = 0; i < stockCloses.size(); i++){
        times.push_back(i);
    }

    double initialValue = stockCloses[0];

    std::vector<std::vector<double>> constants = {{0,0,0}, {0,0,0}};

    //Starting with a very wide constant limit and step values
    std::vector<std::vector<std::vector<double>>> constantLimits = {{{-10,10},{-10,10},{-10,10}},{{-10,10},{-10,10},{-10,10}}};

    std::vector<std::vector<double>> constantSteps = {{0.01,0.01,0.01}, {0.01,0.01,0.01}};

    std::cout << "\nInitial Value: " << initialValue;

    stochasticModel model = stochasticModel(alphaFunction, betaFunction, initialValue, times, constants, constantLimits, constantSteps);

    for(int i = 0; i < 5; i++){

        std::cout << "\nStarting Drift Estimation";
        model.parameters[0] = simulatedAnnealingDriftEstimation(model, 0, stockCloses, 1, 200, 0.9, 500, 1, multiVectorRMSE);
        
        //model.parameters[0] = {0, 0.000865574, 0};
        std::cout << "\nDrift Est: " << model.parameters[0][0] << " " << model.parameters[0][1] << " " << model.parameters[0][2];
        std::cout << "\nStarting Vol Estimation";

        model.betaFunction = betaFunction;
        model.parameters[1] = simulatedAnnealingVolEstimation(model, 1, stockCloses, 500, 50, 0.9, 150, 15, returnComparison);
        
        std::cout << "\n\n For Run " << i;
        std::cout << "\nDrift Est: " << model.parameters[0][0] << " " << model.parameters[0][1] << " " << model.parameters[0][2];;
        std::cout << "\nVolatility Est: " << model.parameters[1][0] << " " << model.parameters[1][1] << " " << model.parameters[1][2];; 

        for(int j = 0; j < constantSteps.size(); j++){
            for(int k = 0; k < constantSteps[j].size(); k++){
                constantSteps[j][k] *= 0.01;
            }
        }

        model.parameterSteps = constantSteps;
    }

    std::vector<std::vector<double>> estimations = multipleEulerMaruyama(model, 500);

    multiVectorToCSV(estimations, "estimationOutput.csv");

    vectorToCSV(stockCloses, "dataOutput.csv");

    return model;
}

polynomialModel fitPolynomial(std::vector<double> &observations, std::vector<double> times, int maxTerm, int divisions){
    
    polynomialModel model = polynomialModel(observations[0], times);

    model.addTerm(1,1);

    std::vector<long double> lowestAIC = {std::numeric_limits<double>::infinity(), 0};

    polynomialModel bestOverall = model;

    while(true){

        int nextTerm = 0;
        std::vector<long double> innerAIC = {std::numeric_limits<double>::infinity(), 0};

        polynomialModel bestModel = model;

        bool notInf = false;

        while(true){
            if(std::find(model.activeTerms[0].begin(), model.activeTerms[0].end(), nextTerm) != model.activeTerms[0].end()){
                nextTerm++;
                continue;
            }
            
            model.addTerm(0,nextTerm);
            model.removeTerm(1,1);
            model.addTerm(1,1);

            for(int i = 0; i < model.parameters[0].size();i++){
                model.parameters[0][i] = 0;
            }
            
            
            model.parameters[0] = polynomialParamEstimation(model, 0, observations, 1, 200, 0.99, 250, 1, driftCost, {double(divisions)});
            
            bool zeroTerm = false;
            for(int i = 0; i < model.activeTerms[0].size(); i++){
                if(model.parameters[0][model.activeTerms[0][i]] == 0.0){
                    zeroTerm = true;
                }
            }

            if(zeroTerm){
                model.removeTerm(0, nextTerm);
                continue;
            }
            

            model.parameters[1] = polynomialParamEstimation(model, 1, observations, 1, 100, 0.9, 100, 1, varianceCost, {double(divisions), 500});

            std::cout << "\nCurrent Model is: \n" << model.toString();

            std::vector<long double> curAIC = model.calculateAIC(observations, 25000, divisions);
            
            std::cout << "\nCurent AIC is: " << curAIC[0] << " , " << curAIC[1];
            //if(!std::isinf(curAIC[0])){
            //    notInf = true;
            //}

            //if(nextTerm < maxTerm && !(notInf)){
            //    bestModel = model;
            //    model.removeTerm(0, nextTerm);
            //    nextTerm++;
            //    continue;
            //}

            //if(curAIC >= innerAIC){
            if(!compareAIC(curAIC, innerAIC)){
                model = bestModel;
                break;
            }else{
                bestModel = model;
                innerAIC = curAIC;
                model.removeTerm(0, nextTerm);
                nextTerm++;
            }
        }

        //if(innerAIC > lowestAIC){
        if(!compareAIC(innerAIC, lowestAIC)){
            model.removeTerm(0, nextTerm);
            break;
        }else{
            lowestAIC = innerAIC;
            bestOverall = model;
        }
    }

    return bestOverall;
}

polynomialModel bestModelNTerms(std::vector<double> &observations, std::vector<double> times, int numTerms, int maxTerm, int divisions){

    polynomialModel model = polynomialModel(polynomialWithXFunction, zeroTimeFunction, observations[0], times);

    model.addMultipleRandomTerms({0,2}, maxTerm, numTerms);

    double startingTemp = 500;
    double tempLimit = 5;
    double coolingRate = 15;
    double stepsAtTemp = 100;
    

    int moveLimit = 1000;

    int notMovedIn = 0;

    double temperature = startingTemp;

    double cost;
    double prob;
    double oldCost = std::numeric_limits<double>::infinity();

    polynomialModel currentModel = model;
    polynomialModel bestModel = model;
    polynomialModel newModel = model;

    double bestCost = std::numeric_limits<double>::infinity();

    while(temperature > tempLimit){

        for(int i = 0; i < stepsAtTemp; i++){
        notMovedIn++;
        
        while(true){
            newModel = polynomialMultiEstimation(newModel, {0,2}, observations, 1, 5000, 0.5, 500, 1, driftCost, {double(divisions), 5000});
            
            bool zeroTerm = false;

            for(int i = 0; i < newModel.activeTerms[0].size(); i++){
                if(newModel.parameters[0][newModel.activeTerms[0][i]] == 0.0){
                    zeroTerm = true;
                }
            }
            for(int i = 0; i < newModel.activeTerms[2].size(); i++){
                if(newModel.parameters[2][newModel.activeTerms[2][i]] == 0.0){
                    zeroTerm = true;
                }
            }
            
            if(!zeroTerm){
                break;
            }
        }
    
        cost = driftCost(newModel, observations, 1, {double(divisions), 5000});


        //std::cout << "\nCur cost: " << cost;
        //std::cout << "\n" << newModel.toString();

        if((cost == 0 || abs(cost) == std::numeric_limits<double>::infinity()) && (oldCost == 0 || abs(oldCost) == std::numeric_limits<double>::infinity())){
            newModel.removeAllTerms();
            newModel.addMultipleRandomTerms({0,2}, maxTerm, numTerms);
            //std::cout << "\nIndeterminate cost!";
            continue;
        }else{
            prob = acceptanceProbability(cost, oldCost, temperature);  
            //std::cout << "\nCur prob: " << prob;     
        }
        
        if(prob > randomGenerator.d01()){
            oldCost = cost;
            currentModel = newModel;
        }

        if(cost < bestCost){
            notMovedIn = 0;
            std::cout << "\nBest cost: " << cost << "\nBest model: \n" << bestModel.toString();
            bestModel = currentModel;
            bestCost = cost;
        }

        if(notMovedIn > moveLimit){
            return bestModel;
        }

        newModel.neighboringSet(maxTerm);

        }
    
        temperature -= coolingRate;
        std::cout << "\nTemp: " << temperature;
    }
    return model;
}

void createPath(){
    auto alphaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * value);
    };

    auto betaFunction = [](double& value, double& time, std::vector<double>& parameters){
        return (parameters[0] * pow(value, parameters[1]));
    };

    double initialValue = 1.0;

    std::vector<std::vector<double>> baseConstants = {{0.05},{0.15,0.5}};

    std::vector<double> times = {0};

    double dt = 0.01;

    while(times.back() < 20){
        times.push_back(times.back() + dt);
    }

    std::cout << times.size() << std::endl;

    stochasticModel obsModel = stochasticModel(alphaFunction, betaFunction, initialValue, times, baseConstants);

    std::vector<std::vector<double>> stockCloses = multipleEulerMaruyama(obsModel, 500);

    multiVectorToCSV(stockCloses, "output.csv");
}