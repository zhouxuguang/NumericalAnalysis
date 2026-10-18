//
// Created by zhouxuguang on 2026/10/17.
// 快速傅里叶变换（FFT）
//

#ifndef NUMERICALANALYSIS_FFT_H
#define NUMERICALANALYSIS_FFT_H

#include <vector>
#include <complex>
#include "PreDefine.h"

using cd = std::complex<double>;

class FFT
{
public:
    FFT(int n);
    
    std::vector<cd> forward(const std::vector<cd>& x) const;

    std::vector<cd> inverse(const std::vector<cd>& X) const;
    
private:
    int mN = 0;
    int mStages = 0;
    std::vector<uint32_t> mReverseTable;
    std::vector<cd> mTw;
    void buildTwiddle();
};

void TestFFT();

#endif //NUMERICALANALYSIS_FFT_H
