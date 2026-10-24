//
// Created by zhouxuguang on 2023/12/2.
//

#include "TestLinearSystem.h"
#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm>

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

    LinearSystem::SORIteration(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}

void TestLinearSystemJacobiEigen()
{
    Matrix A = Matrix(3, 3);
//    A[0][0] = 1;
//    A[0][1] = 3;
//    A[0][2] = 3;
//
//    A[1][0] = -3;
//    A[1][1] = -5;
//    A[1][2] = -3;
//
//    A[2][0] = 3;
//    A[2][1] = 3;
//    A[2][2] = 1;

    A[0][0] = 3;
    A[0][1] = 1;
    A[0][2] = 2;

    A[1][0] = 1;
    A[1][1] = 3;
    A[1][2] = 4;

    A[2][0] = 2;
    A[2][1] = 4;
    A[2][2] = 6;

    Vector X;
    X.resize(3);

    std::vector<Vector> eigenVecotrs;

    EigenValueSystem::JacobiEigenValue(A, X, eigenVecotrs, 0.000000000001, 100);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}

void TestTridiagonalSystem()
{
    Vector A;
    A.push_back(2);
    A.push_back(2);
    A.push_back(-7);
    A.push_back(5);
    //A.push_back(1);

    Vector B;
    B.push_back(1);
    B.push_back(-3);
    B.push_back(4);
    B.push_back(0);
    //B.push_back(0);

    Vector C;
    C.push_back(0);
    C.push_back(1);
    C.push_back(3);
    C.push_back(2);
    //C.push_back(-1);

    Vector F;
    F.push_back(3);
    F.push_back(-3);
    F.push_back(-10);
    F.push_back(2);
    //F.push_back(1);

    Vector X;
    LinearSystem::TridiagonalSystem(A, B, C, F, X);

    printf("追赶法测试结果:\n");
    for (int i = 0; i < X.size(); ++i)
    {
        printf("x[%d] = %lf\n", i, X[i]);
    }
}

void TestLDLT()
{
    Matrix A = Matrix(3, 3);
    A[0][0] = 1;
    A[0][1] = -1;
    A[0][2] = 1;

    A[1][0] = -1;
    A[1][1] = 3;
    A[1][2] = -2;

    A[2][0] = 1;
    A[2][1] = -2;
    A[2][2] = 4.5;

    Vector X;
    X.resize(3);

    Vector B;
    B.push_back(4);
    B.push_back(-8);
    B.push_back(12);

    LinearSystem::LDLTransposeSovle(A, B, X);
    std::cout << "LDLT : x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;
}

void TestLDLT1()
{
    Matrix A = Matrix(4, 4);
    A[0][0] = 4;
    A[0][1] = -2;
    A[0][2] = 4;
    A[0][3] = 2;

    A[1][0] = -2;
    A[1][1] = 10;
    A[1][2] = -2;
    A[1][3] = -7;

    A[2][0] = 4;
    A[2][1] = -2;
    A[2][2] = 8;
    A[2][3] = 4;

    A[3][0] = 2;
    A[3][1] = -7;
    A[3][2] = 4;
    A[3][3] = 7;

    Vector X;
    X.resize(4);

    Vector B;
    B.push_back(8);
    B.push_back(2);
    B.push_back(16);
    B.push_back(6);

    LinearSystem::LDLTransposeSovle(A, B, X);
    std::cout << "LDLT1 : x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << " x4 = " << X[3] << std::endl;
}

//由反射向量v和标量beta构造Householder矩阵 H = I - beta * v * v^T
static void BuildHouseholderMatrix(const Vector &v, double beta, Matrix &H)
{
    int n = (int)v.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            H[i][j] = (i == j ? 1.0 : 0.0) - beta * v[i] * v[j];
        }
    }
}

//验证x的Householder变换：H * x = alpha * e1 且 H 为正交阵
static void CheckHouseholder(const Vector &x)
{
    int n = (int)x.size();

    Vector v;
    double beta = 0.0;
    EigenValueSystem::Householder(x, v, beta);

    Matrix H(n, n);
    BuildHouseholderMatrix(v, beta, H);

    //计算||x||，实现总是把x反射到 +||x|| * e1
    double sigma = 0.0;
    for (int i = 0; i < n; ++i)
    {
        sigma += x[i] * x[i];
    }
    sigma = sqrt(sigma);
    double alpha = sigma;

    printf("Householder: x = (");
    for (int i = 0; i < n; ++i)
    {
        printf("%s%g", (i == 0 ? "" : ", "), x[i]);
    }
    printf("), ||x|| = %g\n", sigma);

    printf("  v = (");
    for (int i = 0; i < n; ++i)
    {
        printf("%s%g", (i == 0 ? "" : ", "), v[i]);
    }
    printf("), beta = %g\n", beta);

    //1、验证 H * x = alpha * e1
    double transformError = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            sum += H[i][j] * x[j];
        }
        double target = (i == 0) ? alpha : 0.0;
        transformError = std::max(transformError, fabs(sum - target));
    }

    //2、验证 H^T * H = I，即H为正交阵
    double orthogonalError = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < n; ++k)
            {
                sum += H[k][i] * H[k][j];
            }
            orthogonalError = std::max(orthogonalError, fabs(sum - (i == j ? 1.0 : 0.0)));
        }
    }

    printf("  max|H*x - alpha*e1| = %e, max|H^T*H - I| = %e\n\n", transformError, orthogonalError);
    assert(transformError < 1e-12);
    assert(orthogonalError < 1e-12);
}

void TestHouseholder()
{
    //2维向量x = (3, 4)^T，H * x 应为 (5, 0)^T
    Vector x2;
    x2.push_back(3);
    x2.push_back(4);
    CheckHouseholder(x2);

    //4维一般向量
    Vector x4;
    x4.push_back(1);
    x4.push_back(-2);
    x4.push_back(3);
    x4.push_back(4);
    CheckHouseholder(x4);

    //尾部元素全为0且x[0] > 0，此时H退化为单位阵
    Vector xPositive;
    xPositive.push_back(5);
    xPositive.push_back(0);
    xPositive.push_back(0);
    CheckHouseholder(xPositive);

    //尾部元素全为0且x[0] < 0，此时H = diag(-1, 1, 1)
    Vector xNegative;
    xNegative.push_back(-5);
    xNegative.push_back(0);
    xNegative.push_back(0);
    CheckHouseholder(xNegative);

    //零向量，此时H退化为单位阵
    Vector x0(3, 0.0);
    CheckHouseholder(x0);
}
