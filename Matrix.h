//
// Created by zhouxuguang on 2023/12/2.
//

#ifndef NUMERICALANALYSIS_MATRIX_H
#define NUMERICALANALYSIS_MATRIX_H

//矩阵类
class Matrix
{
public:
    Matrix(int width, int height);

    ~Matrix();

    double * operator [] (int row);

    /**
     * 交换两行的数据
     * @param row1
     * @param row2
     */
    void SwapRows(int row1, int row2);

    /**
     * 获得第col的元素最大值所在的位置
     * @param col
     * @param startRow
     * @param endRow
     * @return
     */
    int GetColumnMaxAbsElementIndex(int col, int startRow, int endRow);
private:
    int mWidth = 0;
    int mHeight = 0;
    double *mData = nullptr;

};


#endif //NUMERICALANALYSIS_MATRIX_H
