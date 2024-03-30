//
// Created by zhouxuguang on 2024/3/30.
//

#ifndef NUMERICALANALYSIS_INTERPOLATION_H
#define NUMERICALANALYSIS_INTERPOLATION_H

#include <vector>

double GetLagrangeValue(double value, const std::vector<double>& xVec, const std::vector<double>& yVec);

double GetNewtonValue(double value, const std::vector<double>& xVec, const std::vector<double>& yVec);


#endif //NUMERICALANALYSIS_INTERPOLATION_H
