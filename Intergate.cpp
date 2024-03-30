//
//  main.cpp
//  数值积分
//
//  Created by zhouxuguang on 2023/1/1.
//

#include "Intergate.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <math.h>


//复化梯形公式
double ComplexTrapezoidIntergate(double a, double b, int numSection, KernelFunction kernelFunction)
{
    if (a >= b)
    {
        std::swap(a, b);
    }

    double h = (b - a) / numSection;
    double sum = 0.0;
    sum += kernelFunction(a);

    for (int i = 1; i < numSection; i ++)
    {
        sum += (kernelFunction(a + i * h)) * 2;
    }

    sum += kernelFunction(b);
    sum *= h;
    sum *= 0.5;

    return sum;
}

//复化simpson公式求积
double ComplexSimpsonIntergate(double a, double b, int numSection, KernelFunction kernelFunction)
{
    if (a >= b)
    {
        std::swap(a, b);
    }

    numSection *= 2;

    double h = (b - a) / numSection;
    double sum = 0.0;
    sum += kernelFunction(a);

    for (int i = 1; i < numSection; i ++)
    {
        double xk = a + i * h;
        if (i % 2 == 0)   //偶数项
        {
            sum += (kernelFunction(xk)) * 2;
        }
        else
        {
            sum += (kernelFunction(xk)) * 4;
        }
    }

    sum += kernelFunction(b);
    sum *= h;
    sum /= 3;

    return sum;
}

//变步长复化simpson公式求积
double VariableStepComplexSimpsonIntergate(double a, double b, KernelFunction kernelFunction, double eps)
{
    int numSections = 1;
    double s1 = ComplexSimpsonIntergate(a, b, numSections, kernelFunction);

    for (int i = 1; i < 10000; i ++)
    {
        numSections *= 2;
        double s2 = ComplexSimpsonIntergate(a, b, numSections, kernelFunction);
        if (fabs(s2 - s1) < eps * 15)
        {
            return s2;
        }

        s1 = s2;
    }

    return s1;
}

//变步长复化simpson公式求积
double RombergIntergate(double a, double b, KernelFunction kernelFunction, double eps)
{
    const int ttableDims = 10;
    double ttable[ttableDims][ttableDims];  //t值表
    ttable[0][0] = (kernelFunction(a) + kernelFunction(b)) * 0.5 * (b - a);

    int numSections = 1;

    for (int i = 1; i < ttableDims; i ++)
    {
        double h = (b - a) / pow(2.0, i);

        //根据复合梯形公式求积法计算第一列
        ttable[i][0] = ttable[i - 1][0] * 0.5;
        double tempSum = 0;
        for (int k = 1; k <= numSections; k ++)
        {
            tempSum += kernelFunction(a + (2 * k - 1) * h);
        }
        tempSum *= h;

        ttable[i][0] += tempSum;

        //前一次区间二分的次数乘以等于2
        numSections *= 2;

        //计算加速值，即这一行其余的元素值
        for (int k = 1; k <= i; k ++)
        {
            double fourPowK = pow(4, k);
            ttable[i][k] = fourPowK * ttable[i][k - 1] / (fourPowK - 1) - ttable[i - 1][k - 1] / (fourPowK - 1);
        }

        //精度判断
        if (fabs(ttable[i][i] - ttable[i - 1][i - 1]) < eps)
        {
            printf("当前的二分次数 = %d\n", i);
            return ttable[i][i];
        }

    }

    return 0;
}
