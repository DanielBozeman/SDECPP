#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>

void vectorToCSV(std::vector<double> outVector, std::string fileName){

    bool status = std::filesystem::remove(fileName);

    if(status){
        printf( "Successfully removed file\n" );
    } else {
        printf( "No such file!\n"); // No such file or directory
    }

    std::ofstream outfile(fileName);

    if (!outfile.is_open()) {
        std::cout << "Error opening file!\n";
    }else{
        std::cout << "File Opened!\n";
    }

    for(int i = 0; i < outVector.size(); i++){
        outfile << outVector[i] << ",\n";
    }


    std::cout << "File contents written!";

    outfile.close();

}

void multiVectorToCSV(std::vector<std::vector<double>> outVectors, std::string fileName){
    bool status = std::filesystem::remove(fileName);

    if(status){
        printf( "Successfully removed file\n" );
    } else {
        printf( "No such file!\n"); // No such file or directory
    }

    std::ofstream outfile(fileName);

    if (!outfile.is_open()) {
        std::cout << "Error opening file!\n";
    }else{
        std::cout << "File Opened!\n";
    }

    for(int i = 0; i < outVectors[0].size(); i++){
        for(int j = 0; j < outVectors.size(); j++){
            outfile << outVectors[j][i] << ",";
        }
        outfile << std::endl;
    }


    std::cout << "File contents written!";

    outfile.close();
}