//
// Created by zhouxuguang on 2023/12/4.
//

#include "EigenValueSystem.h"
#include <cmath>
#include <iostream>

void
EigenValueSystem::JacobiEigenValue(Matrix &A, Vector &eigenValue, std::vector<Vector> &eigenVectors, double eps,
                                   int maxIterCount)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    int nDim = width;

    int nCount = 0;		//迭代次数
    while (true)
    {
        //在A的非对角线上找到最大元素，以及记录最大元素的位置
        double dbMax = A[0][1];
        int p = 0;   //p
        int q = 1;   //q
        for (int i = 0; i < nDim; i ++)			//行
        {
            for (int j = 0; j < nDim; j ++)		//列
            {
                double d = fabs(A[i][j]);

                if ((i !=j) && (d > dbMax))
                {
                    dbMax = d;
                    p = i;
                    q = j;
                }
            }
        }

        if (dbMax < eps)     //精度符合要求
            break;

        if (nCount > maxIterCount)       //迭代次数超过限制
            break;

        nCount++;

        double dbApp = A[p][p];
        double dbApq = A[p][q];
        double dbAqq = A[q][q];

        double s = (dbAqq - dbApp) / (2 * dbApq);
        double t = 0;
        if (s >= 0.0)
        {
            t = 1.0 / (s + sqrt(1 + s * s));
        }
        else
        {
            t = 1.0 / (s - sqrt(1 + s * s));
        }

        double c = 1.0 / sqrt(1.0 + t * t);
        double d = t * c;

        //计算a[p][q]等相关元素
        A[p][p] = A[p][p] - t * A[p][q];
        A[q][q] = A[q][q] + t * A[p][q];
        A[p][q] = A[q][p] = 0.0;

        //计算第p行和第q行上的元素
        for (int i = 0; i < nDim; i ++)
        {
            if ((i != q) && (i != p))
            {
                double api = A[p][i];
                A[i][p] = A[p][i] = c * A[p][i] - d * A[q][i];  //p行和p列
                A[i][q] = A[q][i] = d * api + c * A[q][i];  //q行和q列
            }
        }

        std::cout << "x1 = " << A[0][0] << " x2 = " <<A[1][1] << " x3 = " << A[2][2] << std::endl;
    }

    //计算特征值
    for (int i = 0; i < nDim; ++i)
    {
        eigenValue[i] = A[i][i];
    }

}
