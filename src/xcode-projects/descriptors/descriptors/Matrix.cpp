//
//  Matrix.cpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#include "Matrix.hpp"

mind::Matrix::Matrix(size_t nrows, size_t ncols, size_t nslices, std::vector<double> elems)
{
    
}

mind::Matrix mind::Matrix::identityMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    
}

mind::Matrix mind::Matrix::onesMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    
}

mind::Matrix mind::Matrix::zerosMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    
}

size_t mind::Matrix::numberOfElements() const
{
    return (numberOfRows_ * numberOfCols_ * numberOfSlices_);
}

double &mind::Matrix::at(size_t row, size_t col, size_t slice)
{
    
}

const double &mind::Matrix::at(size_t row, size_t col, size_t slice) const
{
    
}

size_t mind::Matrix::linearIndexOfElementAtMatrixIndex(size_t row, size_t col, size_t slice) const
{
    return row + col*(numberOfRows_) + slice*(numberOfRows_*numberOfCols_);
}

Matrix mind::Matrix::subtracting(const Matrix &rhs) const
{
    
}

Matrix mind::Matrix::integratedAlongXYZ() const
{
    
}

Matrix mind::Matrix::shiftedBy(const Vector &dr) const
{
    
}

