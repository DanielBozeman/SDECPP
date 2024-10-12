#ifndef __FILEWRITER_H__
#define __FILEWRITER_H__

#include <vector>
#include <string>

void vectorToCSV(std::vector<double> outVector, std::string fileName);

void multiVectorToCSV(std::vector<std::vector<double>> outVectors, std::string fileName);
#endif // __FILEWRITER_H__