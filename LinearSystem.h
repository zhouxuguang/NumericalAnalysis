//
// Created by zhouxuguang on 2023/12/2.
//

#ifndef NUMERICALANALYSIS_LINEARSYSTEM_H
#define NUMERICALANALYSIS_LINEARSYSTEM_H

#include "PreDefine.h"
#include "Matrix.h"

class LinearSystem
{
public:
    /**
     * 列主元高斯消元法
     * @param A
     * @param B
     * @param X
     */
    static void GaussElimination(const Matrix& A, const Vector& B, Vector &X);

    /**
     * 雅可比迭代解线性方程组
     * @param A
     * @param B
     * @param X
     */
    static void JacobiIteration(const Matrix& A, const Vector& B, Vector &X);

    /**
     * 高斯塞德尔迭代求解线性方程组
     * @param A
     * @param B
     * @param X
     */
    static void GaussSeidelIteration(const Matrix& A, const Vector& B, Vector &X);

    /**
     * LDLT分解算法求解对称正定线性方程组
     * @param A
     * @param B
     * @param X
     */
    static void LDLTransposeSovle(const Matrix& A, const Vector& B, Vector &X);

    /**
     * 三对角方程的追赶法，所有的数组元素个数必须相等
     * @param A 主对角线元素
     * @param B 主对角线右侧的元素
     * @param C 主对角线左侧的元素
     * @param f 右端向量
     * @param X 结果数组
     */
    static void TridiagonalSystem(const Vector& A, const Vector &B, const Vector &C, const Vector& f, Vector& X);
};


#endif //NUMERICALANALYSIS_LINEARSYSTEM_H
