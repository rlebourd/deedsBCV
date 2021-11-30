//
//  Matrix.cpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#include "Matrix.hpp"

mind::Matrix::Matrix(size_t nrows, size_t ncols, size_t nslices, std::vector<double> elems)
: numberOfRows_{nrows},
  numberOfCols_{ncols},
  numberOfSlices_{nslices}
{
    assert(elements_.size() == numberOfElements());
    elements_ = elems;
}

mind::Matrix mind::Matrix::fromData(size_t nrows, size_t ncols, size_t nslices, std::vector<double> elements)
{
    return mind::Matrix{nrows, ncols, nslices, elements};
}

mind::Matrix mind::Matrix::identityMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    std::vector<double> elems = std::vector<double>(nrows*ncols*nslices, 0);
    return mind::Matrix{nrows, ncols, nslices, elems};
}

mind::Matrix mind::Matrix::onesMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    std::vector<double> elems = std::vector<double>(nrows*ncols*nslices, 0);
    return mind::Matrix{nrows, ncols, nslices, elems};
}

mind::Matrix mind::Matrix::zerosMatrix(size_t nrows, size_t ncols, size_t nslices)
{
    std::vector<double> elems = std::vector<double>(nrows*ncols*nslices, 0);
    return mind::Matrix{nrows, ncols, nslices, elems};
}

size_t mind::Matrix::numberOfElements() const
{
    return (numberOfRows_ * numberOfCols_ * numberOfSlices_);
}

double &mind::Matrix::at(size_t row, size_t col, size_t slice)
{
    auto linearIndex = linearIndexOfElementAtMatrixIndex(row, col, slice);
    return elements_[linearIndex];
}

const double &mind::Matrix::at(size_t row, size_t col, size_t slice) const
{
    auto linearIndex = linearIndexOfElementAtMatrixIndex(row, col, slice);
    return elements_[linearIndex];
}

double &mind::Matrix::at(const mind::IndexVector &idx)
{
    return at(idx.row(), idx.column(), idx.slice());
}

const double &mind::Matrix::at(const mind::IndexVector &idx) const
{
    return at(idx.row(), idx.column(), idx.slice());
}

size_t mind::Matrix::linearIndexOfElementAtMatrixIndex(size_t row, size_t col, size_t slice) const
{
    assert(row < numberOfRows_);
    assert(col < numberOfCols_);
    assert(slice < numberOfSlices_);

    return row + col*(numberOfRows_) + slice*(numberOfRows_*numberOfCols_);
}

mind::Matrix mind::Matrix::subtracting(const  mind::Matrix &rhs) const
{
    assert(numberOfRows_ == rhs.numberOfRows_);
    assert(numberOfCols_ == rhs.numberOfCols_);
    assert(numberOfSlices_ == rhs.numberOfSlices_);
    
    return identityMatrix(rhs.numberOfRows_, rhs.numberOfCols_, rhs.numberOfSlices_);
}

mind::Matrix mind::Matrix::integratedAlongXYZ() const
{
    return identityMatrix(numberOfRows_, numberOfCols_, numberOfSlices_);
}

mind::Matrix mind::Matrix::translatedBy(const mind::IndexVector &dr) const
{
    return identityMatrix(numberOfRows_, numberOfCols_, numberOfSlices_);
}
