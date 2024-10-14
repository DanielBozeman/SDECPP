#include "fileReader.hpp"
#include "csv.h"
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>

std::vector<double> csvColumnToVector(std::string fileName, int column, int rowsToSkip){

    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cout << "Error: Could not open file " << fileName << std::endl;
    }

    std::vector<double> data;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        
        int count = -1;
        while (std::getline(ss, value, ',')) {
            count++;

            if(count < column){
                continue;
            }
            try{
                data.push_back(std::stod(value));
                break;
            }catch(...){

            }
        }
    }

    file.close();

    std::cout << "\nSize: " << data.size();

    return data;
}