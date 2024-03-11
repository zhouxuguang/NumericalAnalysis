//
// Created by zhouxuguang on 2024/3/10.
//

#ifndef NUMERICALANALYSIS_ODE_H
#define NUMERICALANALYSIS_ODE_H

typedef double (*DiffEquationFunction)(double xValue, double yValue);

class ODE
{
public:
    /**
     * 前向欧拉方法
     * @param function 右端的微分关系式
     * @param x0 x起点
     * @param y0 y起点
     * @param xn x尾点
     * @param n 求解的个数
     */
    static void ForwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n);

    static void BackwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n);

    static void RungeKutta(DiffEquationFunction function, double x0, double y0, double xn, int n);
};


#endif //NUMERICALANALYSIS_ODE_H
