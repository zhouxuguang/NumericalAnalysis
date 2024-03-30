//
// Created by zhouxuguang on 2024/3/30.
//

#ifndef NUMERICALANALYSIS_INTERGATE_H
#define NUMERICALANALYSIS_INTERGATE_H


typedef double (*KernelFunction)(double xValue);

//复化梯形公式
double ComplexTrapezoidIntergate(double a, double b, int numSection, KernelFunction kernelFunction);

//复化simpson公式求积
double ComplexSimpsonIntergate(double a, double b, int numSection, KernelFunction kernelFunction);

//变步长复化simpson公式求积
double VariableStepComplexSimpsonIntergate(double a, double b, KernelFunction kernelFunction, double eps);

//变步长复化simpson公式求积
double RombergIntergate(double a, double b, KernelFunction kernelFunction, double eps);


#endif //NUMERICALANALYSIS_INTERGATE_H
