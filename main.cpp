#include <iostream>
#include "LinearSystem.h"

int main()
{
    Matrix A = Matrix(3, 3);
    A[0][0] = 2;
    A[0][1] = 1;
    A[0][2] = 1;

    A[1][0] = 3;
    A[1][1] = 1;
    A[1][2] = 2;

    A[2][0] = 1;
    A[2][1] = 2;
    A[2][2] = 2;

    Vector X;
    X.resize(3);

    Vector B;
    B.push_back(4);
    B.push_back(6);
    B.push_back(5);

    LinearSystem::GaussElimination(A, B, X);
    std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

    return 0;
}
