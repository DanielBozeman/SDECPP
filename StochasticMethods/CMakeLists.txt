cmake_minimum_required(VERSION 3.5.0)
project(StochasticMethods)

include_directories(../RandomUtils)
if(NOT TARGET RandomUtils)
    add_subdirectory(RandomUtils)
endif()

include_directories(../FileUtils)

include_directories(include)
add_library(${PROJECT_NAME} src/stochasticMethods.cpp src/parameterEstimation.cpp src/modelEstimationTests.cpp)