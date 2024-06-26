//
// Created by zhouxuguang on 2024/3/10.
//

#include "ODE.h"
#include <cstdio>
#include <cmath>

void ODE::ForwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n)
{
    double h = (xn - x0) / n;
    double x = x0;
    double y = y0;

    int i = 0;
    printf("向前欧拉方法\n");
    printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    for (i = 1;  i <= n; i ++)
    {
        y = y + h * function(x, y);
        x += h;
        printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    }
}

void ODE::BackwardEuler(DiffEquationFunction function, double x0, double y0, double xn, int n)
{
    double h = (xn - x0) / n;
    double x = x0;
    double y = y0;

    int i = 0;
    printf("向后欧拉方法\n");
    printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    for (i = 1;  i <= n; i ++)
    {
        double prevY = y;
        double tempY = y + h * function(x, y);

        //更新之后的x值
        x += h;

        //隐式格式的迭代
        while (true)
        {
            y = prevY + h * function(x, tempY);

            //前后两次迭代的值十分接近，停止迭代
            if (abs(tempY - y) <= 1E-10)
            {
                break;
            }
            tempY = y;
        }

        printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    }
}

void ODE::RungeKutta(DiffEquationFunction function, double x0, double y0, double xn, int n)
{
    double h = (xn - x0) / n;
    double x = x0;
    double y = y0;

    int i = 0;
    printf("Runge-Kutta方法\n");
    printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    for (i = 1;  i <= n; i ++)
    {
        double k1 = function(x, y);
        double k2 = function(x + 0.5 * h, y + 0.5 * h * k1);
        double k3 = function(x + 0.5 * h, y + 0.5 * h * k2);
        double k4 = function(x + h, y + h * k3);

        y = y + (h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0);
        x += h;
        printf("i = %d, x = %lf, y = %lf\n", i, x, y);
    }
}
