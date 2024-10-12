cmake_minimum_required(VERSION 3.5.0)
project(RandomUtils)

include_directories(include)
add_library(${PROJECT_NAME} src/brownianPaths.cpp RandomUtils.hpp)