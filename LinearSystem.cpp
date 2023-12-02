//
// Created by zhouxuguang on 2023/12/2.
//

#include "LinearSystem.h"

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

    //每一行逐步进行消元
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
