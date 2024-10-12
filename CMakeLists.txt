cmake_minimum_required(VERSION 3.5.0)
project(MainProject)

include_directories(RandomUtils)
if(NOT TARGET RandomUtils)
    add_subdirectory(RandomUtils)
endif()

include_directories(StochasticMethods)
if(NOT TARGET StochasticMethods)
    add_subdirectory(StochasticMethods)
endif()

include_directories(FileUtils)
if(NOT TARGET FileUtils)
    add_subdirectory(FileUtils)
endif()

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} RandomUtils StochasticMethods FileUtils)