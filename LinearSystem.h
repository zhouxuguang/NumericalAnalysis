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
     * 选主元高斯消元法
     * @param A
     * @param B
     * @param X
     */
    static void GaussElimination(const Matrix& A, const Vector& B, Vector &X);
};


#endif //NUMERICALANALYSIS_LINEARSYSTEM_H
