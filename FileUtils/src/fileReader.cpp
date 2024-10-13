#include "fileReader.hpp"
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>

void csvColumnToVector(std::string fileName, int column){

    std::ifstream infile(fileName);

    if (!infile.is_open()) {
        std::cout << "\nError opening file!";
    }else{
        std::cout << "\nFile Opened!";
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::cout << line << std::endl;
    }

    std::cout << "\nFile contents read!";

    infile.close();

    return;
}