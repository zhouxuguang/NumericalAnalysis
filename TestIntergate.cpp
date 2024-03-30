//
// Created by zhouxuguang on 2024/3/30.
//

#include "TestIntergate.h"
#include "Intergate.h"

#include <stdio.h>

double pai_intergate(double value)
{
    return 1.0 / (1.0 + value * value);
}

void TestIntergate()
{
    double sum1 = ComplexTrapezoidIntergate(0, 1, 128, pai_intergate);
    sum1 *= 4;

    double sum2 = ComplexSimpsonIntergate(0, 1, 8, pai_intergate);
    sum2 *= 4;

    double sum3 = VariableStepComplexSimpsonIntergate(0, 1, pai_intergate, 1E-7 * 0.5);
    sum3 *= 4;

    double sum4 = RombergIntergate(0, 1, pai_intergate, 1E-7 * 0.5);
    sum4 *= 4;

    printf("复化梯形公式 = %lf, 复化simpson公式 = %lf, 变步长simpson公式 = %lf, Romberg公式 = %lf\n", sum1, sum2, sum3, sum4);
}
