#include "brownianPaths.hpp"
#include <math.h>
#include <chrono>

double inverseNormal(double p){
double a1 = -39.69683028665376;
double a2 =  220.9460984245205;
double a3 = -275.9285104469687;
double a4 =  138.3577518672690;
double a5 = -30.66479806614716;
double a6 =  2.506628277459239;

double b1 = -54.47609879822406;
double b2 =  161.5858368580409;
double b3 = -155.6989798598866;
double b4 =  66.80131188771972;
double b5 = -13.28068155288572;

double c1 = -0.007784894002430293;
double c2 = -0.3223964580411365;
double c3 = -2.400758277161838;
double c4 = -2.549732539343734;
double c5 =  4.374664141464968;
double c6 =  2.938163982698783;

double d1 =  0.007784695709041462;
double d2 =  0.3224671290700398;
double d3 =  2.445134137142996;
double d4 =  3.754408661907416;

double p_low  = 0.02425;
double p_high = 1 - p_low;

double x = 0;

if(0 < p && p < p_low){
   double q = sqrt(-2*log(p));
   x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
}

if(p_low <= p && p <= p_high){
   double q = p - 0.5;
   double r = q*q;
   x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
}

if (p_high < p && p < 1){
   double q = sqrt(-2*log(1-p));
   x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
}

return x;
}
Xoshiro256plus timeSeedRand()
{
    return(Xoshiro256plus(int(std::chrono::system_clock::now().time_since_epoch().count())));
} 

Xoshiro256plus randomGenerator = timeSeedRand();

randomPathMaker::randomPathMaker(){
}

double randomPathMaker::dW(double dt){
    return (inverseNormal(randomGenerator.d01())* dt);
}

std::vector<double> randomPathMaker::makePath(double intervalStart, double intervalEnd, double timeStep){
            
            std::vector<double> brownianPath = {0};

            brownianPath.reserve((intervalEnd - intervalStart)/timeStep);

            double totalInterval = intervalEnd - intervalStart;

            double currentTime = intervalStart;

            while (currentTime < intervalEnd-timeStep){
                brownianPath.push_back((brownianPath.back() + dW(timeStep)));
                currentTime += timeStep;
            }

            return brownianPath;
        }

std::vector<std::vector<double>> randomPathMaker::makeMultiplePaths(double intervalStart, double intervalEnd, double timeStep, int numPaths){
            std::vector<std::vector<double>> paths;

            paths.reserve(numPaths * ((intervalEnd - intervalStart)/timeStep));
            
            for(int i = 0; i < numPaths; i++){
                std::vector<double> brownianPath = {0};

                brownianPath.reserve((intervalEnd - intervalStart)/timeStep);

                double totalInterval = intervalEnd - intervalStart;

                double currentTime = intervalStart;

                while (currentTime < intervalEnd-timeStep){
                    brownianPath.push_back((brownianPath.back() + dW(timeStep)));
                    currentTime += timeStep;
                }

                paths.push_back(brownianPath);
            }

            return paths;
        }

std::vector<std::vector<double>> randomPathMaker::makeCorrelatedPaths(double intervalStart, double intervalEnd, double timeStep, double correlation){
    std::vector<double> path1 = makePath(intervalStart, intervalEnd, timeStep);
    std::vector<double> path2 = makePath(intervalStart, intervalEnd, timeStep);

    std::vector<double> correlatedPath;

    for(int i = 0; i < path1.size(); i++){
        correlatedPath.push_back( correlation * path1[i] + sqrt(1 - correlation*correlation ) * path2[i]);
    }

    return {path1, correlatedPath};
}




