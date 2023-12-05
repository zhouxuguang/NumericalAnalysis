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
};


#endif //NUMERICALANALYSIS_EIGENVALUESYSTEM_H
