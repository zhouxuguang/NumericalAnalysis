//
// Created by zhouxuguang on 2026/10/17.
// 离散傅里叶变换（DFT）
//

#include "DFT.h"
#include <algorithm>
#include <math.h>

DFT::DFT(int n) : mN(n), mTw(n)
{
    // 构造旋转矩阵，根据旋转矩阵的周期性，这里只需要存储n个数值就可以了
    const double PI = std::acos(-1.0);
    
    for (int i = 0; i < mN; i ++)
    {
        double theta = -2.0 * PI * i / mN;
        mTw[i] = cd(cos(theta), sin(theta));
    }
}

std::vector<cd> DFT::forward(const std::vector<cd>& x) const
{
    std::vector<cd> result(mN, cd(0.0, 0.0));
    for (int k = 0; k < mN; k ++)
    {
        for (int n = 0; n < mN; n ++)
        {
            const long long idx = static_cast<long long>(k) * n % mN;
            result[k] += mTw[idx] * x[n];
        }
    }
    
    return result;
}

std::vector<cd> DFT::inverse(const std::vector<cd>& X) const
{
    std::vector<cd> result(mN, cd(0.0, 0.0));
    for (int k = 0; k < mN; k ++)
    {
        for (int n = 0; n < mN; n ++)
        {
            const long long idx = static_cast<long long>(k) * n % mN;
            result[k] += std::conj(mTw[idx]) * X[n];
        }
        result[k] /= mN;
    }
    
    return result;
}

void TestDFT()
{
    std::vector<cd> input;
    input.push_back(cd(2.3));
    input.push_back(cd(3.0));
    input.push_back(cd(3.8));
    input.push_back(cd(4.2));
    input.push_back(cd(5.1));
    input.push_back(cd(4.5));
    input.push_back(cd(3.2));
    input.push_back(cd(1.8));
    
    DFT dft(input.size());
    std::vector<cd> output = dft.forward(input);
    
    printf("dft forward : \n");
    for (const auto& iter : output)
    {
        printf("real = %lf, imag = %lf\n", iter.real(), iter.imag());
    }
    
    std::vector<cd> inverseOutput = dft.inverse(output);
    
    printf("dft inverse : \n");
    for (const auto& iter : inverseOutput)
    {
        printf("real = %lf, imag = %lf\n", iter.real(), iter.imag());
    }
}
