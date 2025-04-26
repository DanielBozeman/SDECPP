#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>

void vectorToCSV(std::vector<double> outVector, std::string fileName){

    bool status = std::filesystem::remove(fileName);

    if(status){
        printf( "\nSuccessfully removed file" );
    } else {
        printf( "\nNo such file!"); // No such file or directory
    }

    std::ofstream outfile(fileName);

    if (!outfile.is_open()) {
        std::cout << "\nError opening file!";
    }else{
        std::cout << "\nFile Opened!";
    }

    for(int i = 0; i < outVector.size(); i++){
        outfile << outVector[i] << ",\n";
    }


    std::cout << "\nFile contents written!";

    outfile.close();

}

void multiVectorToCSV(std::vector<std::vector<double>> outVectors, std::string fileName){
    bool status = std::filesystem::remove(fileName);

    if(status){
        //printf( "\nSuccessfully removed file" );
    } else {
        //printf( "\nNo such file!"); // No such file or directory
    }

    std::ofstream outfile(fileName);

    if (!outfile.is_open()) {
        std::cout << "\nError opening file!";
    }else{
        //std::cout << "\nFile Opened!";
    }

    for(int i = 0; i < outVectors[0].size(); i++){
        for(int j = 0; j < outVectors.size(); j++){
            outfile << outVectors[j][i] << ",";
        }
        outfile << std::endl;
    }


    std::cout << "\nFile contents written!";

    outfile.close();
}