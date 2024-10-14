#ifndef __FILEREADER_H__
#define __FILEREADER_H__

#include <string>
#include <vector>

std::vector<double> csvColumnToVector(std::string fileName, int column, int rowsToSkip = 0);

#endif // __FILEREADER_H__