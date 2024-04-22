cmake_minimum_required(VERSION 3.22)
project(BPO-SMK)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(.)
include_directories(SFMT)
include_directories(IM)
include_directories(Graph)
include_directories(RISCalculation)
include_directories(ParetoOptimization)


file(GLOB Graph_SOURCES
        "Graph/*.hpp"
        "Graph/*.cpp")

file(GLOB RISCalculation_SOURCES
        "RISCalculation/*.hpp"
        "RISCalculation/*.cpp")

file(GLOB IM_SOURCES
        "IM/*.hpp"
        "IM/*.cpp"
        "BasicGreedy/*.hpp"
        "BasicGreedy/*.cpp")

file(GLOB SFMT_SOURCES
        "SFMT/SFMT.h"
        "SFMT/SFMT.c")

set(Common_SOURCES
        ${SFMT_SOURCES}
        log.h
        memoryusage.h
        json.hpp
        IM/TIMInfluenceCalculator.hpp
        IM/TIMInfluenceCalculator.cpp
        IM/TIMCoverage.hpp
        IM/TIMCoverage.cpp
        IM/TIMNormalizedCoverage.cpp
        IM/TIMNormalizedCoverage.hpp
        ParetoOptimization/BPO.h
        ParetoOptimization/BPO.cpp
        ParetoOptimization/ModularCostGenerator.cpp
        ParetoOptimization/ModularCostGenerator.hpp
        RISCalculation/InfluenceCalcSpaceEfficient.cpp
        RISCalculation/InfluenceCalcSpaceEfficient.hpp
        RISCalculation/RISCalcSpaceEfficient.hpp
        RISCalculation/RISCalcSpaceEfficient.cpp

        Graph/Graph.cpp
        Graph/Graph.hpp
        ParetoOptimization/BloomFilter.h
        ParetoOptimization/BloomFilter.cpp
        SensorPlacement/Sensor.h
        SensorPlacement/Sensor.cpp
        ParetoOptimization/DynamicBPO.h
        ParetoOptimization/DynamicBPO.cpp)

add_executable(IM
        ${Common_SOURCES}
        ${IM_SOURCES}
        ${Graph_SOURCES}
        main.cpp)


