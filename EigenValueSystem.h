//
// Created by zhouxuguang on 2023/12/4.
//

#ifndef NUMERICALANALYSIS_EIGENVALUESYSTEM_H
#define NUMERICALANALYSIS_EIGENVALUESYSTEM_H

#include "PreDefine.h"
#include "Matrix.h"

class EigenValueSystem
{
public:
    /**
     * 实对称矩阵特征值和特征向量的雅可比迭代算法
     * @param A
     * @param eigenValue
     * @param eigenVectors
     * @param eps
     * @param maxIterCount
     */
    static void JacobiEigenValue(const Matrix& A,
                                 Vector & eigenValue,
                                 std::vector<Vector>& eigenVectors,
                                 double eps,
                                 int maxIterCount);

    /**
     * 计算向量x的Householder变换（初等反射）
     * 求单位化前的反射向量v和标量beta，使得
     * H = I - beta * v * v^T 满足 H * x = ||x|| * e1
     * 注意：若x为零向量，则返回 v = 0, beta = 0，此时H为单位阵
     * @param x     输入向量（长度n，且n > 0）
     * @param v     输出的Householder向量（长度n）
     * @param beta  输出标量，H = I - beta * v * v^T
     */
    static void Householder(const Vector& x, Vector& v, double& beta);
};


#endif //NUMERICALANALYSIS_EIGENVALUESYSTEM_H
