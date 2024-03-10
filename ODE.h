//
// Created by zhouxuguang on 2024/3/10.
//

#ifndef NUMERICALANALYSIS_ODE_H
#define NUMERICALANALYSIS_ODE_H

typedef double (*DiffEquationFunction)(double xValue, double yValue);

class ODE
{
public:
    static void ForwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n);

    static void BackwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n);
};


#endif //NUMERICALANALYSIS_ODE_H
