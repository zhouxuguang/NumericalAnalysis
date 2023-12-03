//
// Created by zhouxuguang on 2023/12/2.
//

#include "TestLinearSystem.h"

void TestLinearSystem1()
{
    Matrix A = Matrix(3, 3);
    A[0][0] = 2;
    A[0][1] = 1;
    A[0][2] = 1;

    A[1][0] = 3;
    A[1][1] = 1;
    A[1][2] = 2;

    A[2][0] = 1;
    A[2][1] = 2;
    A[2][2] = 2;

    Vector X;
    X.resize(3);

    Vector B;
    B.push_back(4);
    B.push_back(6);
    B.push_back(5);

    LinearSystem::GaussElimination(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}

//这个测试数值稳定性
void TestLinearSystem2()
{
    Matrix A = Matrix(2, 2);
    A[0][0] = 1e-19;
    A[0][1] = 1;

    A[1][0] = 1;
    A[1][1] = 1;


    Vector X;
    X.resize(2);

    Vector B;
    B.push_back(1);
    B.push_back(2);

    LinearSystem::GaussElimination(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << std::endl;
}

void TestLinearSystemJacobi()
{
    Matrix A = Matrix(3, 3);
    A[0][0] = 10;
    A[0][1] = -1;
    A[0][2] = -2;

    A[1][0] = -1;
    A[1][1] = 10;
    A[1][2] = -2;

    A[2][0] = -1;
    A[2][1] = -1;
    A[2][2] = 5;

    Vector X;
    X.resize(3);

    Vector B;
    B.push_back(7.2);
    B.push_back(8.3);
    B.push_back(4.2);

    LinearSystem::JacobiIteration(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}

void TestLinearSystemGaussSidel()
{
    Matrix A = Matrix(3, 3);
    A[0][0] = 10;
    A[0][1] = -1;
    A[0][2] = -2;

    A[1][0] = -1;
    A[1][1] = 10;
    A[1][2] = -2;

    A[2][0] = -1;
    A[2][1] = -1;
    A[2][2] = 5;

    Vector X;
    X.resize(3);

    Vector B;
    B.push_back(7.2);
    B.push_back(8.3);
    B.push_back(4.2);

    LinearSystem::GaussSidelIteration(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}
