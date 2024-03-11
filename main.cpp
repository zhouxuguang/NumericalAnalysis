#include <iostream>
#include "LinearSystem.h"
#include "TestLinearSystem.h"
#include "ODE.h"

int clz(uint32_t x)
{
    int n = 0;
    if (x <= 0x0000ffff)
    {
        n += 16;
        x <<= 16;
    }
    if (x <= 0x00ffffff)
    {
        n += 8;
        x <<= 8;
    }
    if (x <= 0x0fffffff)
    {
        n += 4;
        x <<= 4;
    }
    if (x <= 0x3fffffff)
    {
        n += 2;
        x <<= 2;
    }
    if (x <= 0x7fffffff)
    {
        n ++;
    }

    return n;
}

double f1(double x, double y)
{
    return x - y + 1;
}

double f2(double x, double y)
{
    return x * x * x - (y / x);
}

int main()
{
    //TestLinearSystemGaussSidel();
    TestLinearSystemJacobiEigen();

    uint32_t a = 20;
    uint32_t base = 31 - clz(a);
    printf("log : %d = %d\n", a, base);

    TestTridiagonalSystem();
    TestLDLT();
    TestLDLT1();

    ODE::ForwardEuler(f1, 0, 1, 0.5, 5);
    ODE::BackwardEuler(f1, 0, 1, 0.5, 5);
    ODE::RungeKutta(f2, 1.0, 0.4, 1.9, 9);

    return 0;
}
