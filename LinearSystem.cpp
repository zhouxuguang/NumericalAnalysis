//
// Created by zhouxuguang on 2023/12/2.
//

#include "LinearSystem.h"
#include <limits>
#include <iostream>
#include <assert.h>

void LinearSystem::GaussElimination(const Matrix &A, const Vector &B, Vector &X)
{
    //1 构造增广矩阵
    Matrix extendMat = Matrix(A.GetWith() + 1, A.GetHeight());
    for (int i = 0; i < A.GetHeight(); ++i)
    {
        memcpy(extendMat[i], A[i], A.GetWith() * sizeof(double));
        extendMat[i][A.GetWith()] = B[i];
    }

    int height = A.GetHeight();
    int width = extendMat.GetWith();

    //每一行逐步进行消元，即遍历对角元
    for (int k = 0; k < height - 1; ++k)
    {
        int maxIndex = extendMat.GetColumnMaxAbsElementIndex(k, k, height - 1);

        //当前行的第一个元素，akk不是绝对值最大元素，就需要将当前元素所在行和绝对值最大的行进行交换
        if (maxIndex != k)
        {
            extendMat.SwapRows(k, maxIndex);
        }

        //对第k行下方的所有行进行遍历
        for (int i = k + 1; i < height; ++i)
        {
            //计算好每一行的需要乘以的系数
            extendMat[i][k] = extendMat[i][k] / extendMat[k][k];

            //对第i行，从k+1列开始遍历
            for (int j = k + 1; j < width; ++j)
            {
                extendMat[i][j] = extendMat[i][j] - extendMat[i][k] * extendMat[k][j];
            }
        }
    }

    //开始回代求解
    for (int i = height - 1; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j < height; ++j)
        {
            sum += extendMat[i][j] * X[j];
        }

        X[i] = (extendMat[i][width - 1] - sum) / extendMat[i][i];
    }
}

void LinearSystem::JacobiIteration(const Matrix &A, const Vector &B, Vector &X)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    Vector X0;
    X0.resize(X.size());
    memset(X0.data(), 0, X0.size() * sizeof(double));

    Vector diffVec;
    diffVec.resize(X.size());

    while (true)
    {
        for (int i = 0; i < width; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < width; ++j)
            {
                if (i == j)
                {
                    continue;
                }
                sum += A[i][j] * X0[j];
            }

            // xi(k1)a ( ax(k) ax(k)b)
            X[i] = (B[i] - sum) / A[i][i];
        }

        //判断X0和X相差的阈值小于一定量，停止迭代
        double sum = 0.0;
        for (int i = 0; i < width; ++i)
        {
            diffVec[i] = X[i] - X0[i];
            sum += diffVec[i] * diffVec[i];
        }

        std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

        if (sum <= std::numeric_limits<double>::epsilon())
        {
            break;
        }

        //更新X0
        memcpy(X0.data(), X.data(), X.size() * sizeof(double));
    }
}

void LinearSystem::GaussSeidelIteration(const Matrix &A, const Vector &B, Vector &X)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    Vector X0;
    X0.resize(X.size());
    memset(X0.data(), 0, X0.size() * sizeof(double));

    Vector diffVec;
    diffVec.resize(X.size());

    while (true)
    {
        for (int i = 0; i < width; ++i)
        {
            double sum = 0.0;

            //在计算xi的时候，新的xi前面的都计算好了，直接用新的来计算
            for (int j = 0; j < i; ++j)
            {
                sum += A[i][j] * X[j];
            }

            for (int j = i + 1; j < width; ++j)
            {
                sum += A[i][j] * X0[j];
            }

            X[i] = (B[i] - sum) / A[i][i];
        }

        //判断X0和X相差的阈值小于一定量，停止迭代
        double sum = 0.0;
        for (int i = 0; i < width; ++i)
        {
            diffVec[i] = X[i] - X0[i];
            sum += diffVec[i] * diffVec[i];
        }

        std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

        if (sum <= std::numeric_limits<double>::epsilon())
        {
            break;
        }

        //更新X0
        memcpy(X0.data(), X.data(), X.size() * sizeof(double));
    }
}

void LinearSystem::SORIteration(const Matrix& A, const Vector& B, Vector &X)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    // x0 代表上一次的迭代结果
    Vector X0;
    X0.resize(X.size());
    memset(X0.data(), 0, X0.size() * sizeof(double));

    Vector diffVec;
    diffVec.resize(X.size());
    
    int iterCount = 1;
    double worOpt = 1.1;

    while (true)
    {
        for (int i = 0; i < width; ++i)
        {
            double sum = 0.0;

            //在计算xi的时候，新的xi前面的都计算好了，直接用新的来计算
            for (int j = 0; j < i; ++j)
            {
                sum += A[i][j] * X[j];
            }

            for (int j = i + 1; j < width; ++j)
            {
                sum += A[i][j] * X0[j];
            }

            X[i] = X0[i] * (1.0 - worOpt) + worOpt * (B[i] - sum) / A[i][i];
        }

        //判断X0和X相差的阈值小于一定量，停止迭代
        double sum = 0.0;
        for (int i = 0; i < width; ++i)
        {
            diffVec[i] = X[i] - X0[i];
            sum += diffVec[i] * diffVec[i];
        }

        std::cout <<"iterCount = " << iterCount << ": x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

        if (sum <= std::numeric_limits<double>::epsilon())
        {
            break;
        }

        //更新X0
        memcpy(X0.data(), X.data(), X.size() * sizeof(double));
        iterCount ++;
    }
}

void LinearSystem::LDLTransposeSovle(const Matrix &A, const Vector &B, Vector &X)
{
    assert(A.GetWith() == A.GetHeight());
    int n = A.GetWith();
    Matrix L(n, n);

    Vector D(n);
    Vector Y(n);
    Vector Z(n);
    X.resize(n);

    //分解步骤
    for (int k = 0; k < n; ++k)
    {
        double sum = 0.0;
        for (int r = 0; r < k; ++r)
        {
            sum += D[r] * L[k][r] * L[k][r];
        }
        D[k] = A[k][k] - sum;
        //printf("D[%d] = %lf\n", k + 1, D[k]);

        for (int i = k + 1; i < n; ++i)
        {
            double sum = 0.0;
            for (int r = 0; r < k; ++r)
            {
                sum += D[r] * L[i][r] * L[k][r];
            }

            //printf("sum = %lf\n", sum);

            L[i][k] = (A[i][k] - sum) / D[k];
            //printf("L[%d][%d] = %lf\n", i + 1, k + 1, L[i][k]);
        }
    }

    //回代步骤
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
        {
            sum += L[i][j] * Z[j];
        }
        Z[i] = B[i] - sum;
        Y[i] = Z[i] / D[i];
    }

    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += L[j][i] * X[j];
        }

        X[i] = Y[i] - sum;
    }
}

void LinearSystem::TridiagonalSystem(const Vector &A, const Vector &B, const Vector &C, const Vector &f, Vector &X)
{
    assert(A.size() == B.size());
    assert(A.size() == C.size());
    assert(A.size() == f.size());

    size_t count = A.size();

    //初始化
    Vector u(count);
    Vector v(count);
    Vector y(count);

    u[0] = A[0];
    v[0] = B[0] / u[0];
    y[0] = f[0] / u[0];

    //追的循环
    for (int i = 1; i < count; ++i)
    {
        u[i] = A[i] - C[i] * v[i - 1];
        v[i] = B[i] / u[i];
        y[i] = (f[i] - C[i] * y[i - 1]) / u[i];
    }

    X.resize(count);
    X[count - 1] = y[count - 1];

    //赶的循环
    for (int i = (int)count - 2; i >= 0; --i)
    {
        X[i] = y[i] - v[i] * X[i+1];
    }
}
