//
// Created by zhouxuguang on 2024/3/30.
//

#include "Interpolation.h"

//
//  main.cpp
//  Lagrange插值
//
//  Created by zhouxuguang on 2022/12/12.
//

#include <iostream>
#include <vector>

double GetLagrangeValue(double value, const std::vector<double>& xVec, const std::vector<double>& yVec)
{
    //检查数据个数
    if (xVec.size() != yVec.size())
    {
        return 0.0;
    }

    double sum = 0.0;

    for (int i = 0; i < xVec.size(); i ++)
    {
        double up = 1;   //分子
        double low = 1;   //分母
        for (int j = 0; j < xVec.size(); j ++)
        {
            if (j != i)
            {
                up *= (value - xVec[j]);
                low *= (xVec[i] - xVec[j]);
            }
        }

        sum += yVec[i] * (up / low);
    }

    return sum;
}

//https://www.bilibili.com/read/cv16912961?spm_id_from=333.999.0.0
double GetNewtonValue(double value, const std::vector<double>& xVec, const std::vector<double>& yVec)
{
    //检查数据个数
    if (xVec.size() != yVec.size())
    {
        return 0.0;
    }

    //构造牛顿差商表
    double *diff = new double[xVec.size() * xVec.size()];
    memcpy(diff, yVec.data(), yVec.size() * sizeof(double));
    for (int i = 1; i < xVec.size(); i ++)
    {
        for (int j = 0; j < xVec.size(); j ++)
        {
            diff[i * xVec.size() + j] = (diff[(i - 1) * xVec.size() + j + 1] - diff[(i - 1) * xVec.size() + j]) /
                                        (xVec[i + j] - xVec[j]);
        }
    }

    //最后使用查商标的第一列就可以了
    double sum = 0.0;
    double polyNum = 1.0;  //多项式项的值
    for (int i = 0; i < xVec.size(); i ++)
    {
        sum += polyNum * diff[i * xVec.size()];
        polyNum *= (value - xVec[i]);         //累乘得到下一次多项式的项
    }

    delete []diff;

    return sum;
}
