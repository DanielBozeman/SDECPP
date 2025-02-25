cmake_minimum_required(VERSION 3.10)
project(FileUtils)

include_directories(include)
add_library(${PROJECT_NAME} src/fileWriter.cpp src/fileReader.cpp)
