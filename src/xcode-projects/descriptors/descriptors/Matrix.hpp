//
//  Matrix.hpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>

namespace mind {

    class Matrix {
        const size_t numberOfRows_;
        const size_t numberOfCols_;
        const size_t numberOfSlices_;
        const std::vector<double> elements_;
        
        Matrix(size_t nrows, size_t ncols, size_t nslices, std::vector<double> elems);
        
        size_t linearIndexOfElementAtMatrixIndex(size_t row, size_t col, size_t slice) const;
        
    public:
        static Matrix identityMatrix(size_t nrows, size_t ncols, size_t nslices);
        static Matrix onesMatrix(size_t nrows, size_t ncols, size_t nslices);
        static Matrix zerosMatrix(size_t nrows, size_t ncols, size_t nslices);

        size_t numberOfElements() const;
        double &at(size_t row, size_t col, size_t slice);
        const double &at(size_t row, size_t col, size_t slice) const;
        Matrix subtracting(const Matrix &rhs) const;
        Matrix integratedAlongXYZ() const;
        Matrix shiftedBy(const Vector &dr) const;
        
        /// @brief A Callable object is a function-like object that accepts
        /// a matrix index and element value as its arguments, and which
        /// returns the transformed element
        ///
        /// The resulting matrix is constructed by applying the given
        /// transformation to each element in the matrix.
        template <typename Callable>
        Matrix mappingElements(Callable &&func) const;
    };

}

#endif /* Matrix_hpp */
