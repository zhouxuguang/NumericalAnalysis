//
// Created by zhouxuguang on 2026/10/17.
// 离散傅里叶变换（DFT）
//

#ifndef NUMERICALANALYSIS_DFT_H
#define NUMERICALANALYSIS_DFT_H

#include <vector>
#include <complex>
#include "PreDefine.h"

using cd = std::complex<double>;

class DFT
{
public:
    DFT(int n);
    
    // 正变换
    std::vector<cd> forward(const std::vector<cd>& x) const;
    
    // 逆变换
    std::vector<cd> inverse(const std::vector<cd>& X) const;
    
private:
    int mN;
    std::vector<cd> mTw;
};

void TestDFT();

#endif //NUMERICALANALYSIS_DFT_H
