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
    std::vector<double> brown;
    randomPathMaker rp;

    std::string fileName = "output.csv";

    int status = remove("output.csv");

    std::ofstream outfile(fileName);

    auto t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < 1; i++){
        std::vector<std::vector<double>> paths = rp.makeCorrelatedPaths(0, 10, 0.005, -0.5);

        std::vector<double> brown1 = paths[0];
        std::vector<double> brown2 = paths[1];

        for(int i = 0; i < paths[0].size(); i++){
            for(int j = 0; j < paths.size(); j++){
                outfile << paths[j][i] << ",";
            }
            outfile << "\n";
        }
    }
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

std::vector<std::vector<double>> eulerMaruyamaTest(int numSims){

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

    std::vector<std::vector<double>> approximations;

    for( int i = 0; i < numSims; i++){
        approximations.push_back(eulerMaruyama(alphaFunction, betaFunction, initialValue, times, constants));
    }
    return approximations;
}

void fileWriterTest(){

    std::vector<std::vector<double>> vectors = {};

    auto t1 = std::chrono::high_resolution_clock::now();

    
    eulerMaruyamaTest(10000);

    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "\nTime taken: " << ms_double.count()/1000;
}

int main(){
    fileWriterTest();
}