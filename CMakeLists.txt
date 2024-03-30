cmake_minimum_required(VERSION 3.26)
project(NumericalAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(NumericalAnalysis main.cpp
        Matrix.cpp
        Matrix.h
        PreDefine.h
        LinearSystem.cpp
        LinearSystem.h
        TestLinearSystem.cpp
        TestLinearSystem.h
        EigenValueSystem.cpp
        EigenValueSystem.h
        ODE.cpp
        ODE.h
        Interpolation.cpp
        Interpolation.h
        Intergate.cpp
        Intergate.h
        TestIntergate.cpp
        TestIntergate.h)
