cmake_minimum_required(VERSION 3.26)
project(NumericalAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(NumericalAnalysis main.cpp
        Matrix.cpp
        Matrix.h
        PreDefine.h
        LinearSystem.cpp
        LinearSystem.h)
