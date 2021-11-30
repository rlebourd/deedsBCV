//
//  IndexVector.cpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#include "IndexVector.hpp"

mind::IndexVector::IndexVector(size_t row, size_t col, size_t slice)
: row_{row},
  col_{col},
  slice_{slice}
{
    
}

mind::IndexVector mind::IndexVector::subtracting(const IndexVector &rhs) const
{
    return IndexVector{row_ - rhs.row_, col_ - rhs.col_, slice_ - rhs.slice_};
}

mind::IndexVector mind::IndexVector::adding(const IndexVector &rhs) const
{
    return IndexVector{row_ + rhs.row_, col_ + rhs.col_, slice_ + rhs.slice_};
}

size_t mind::IndexVector::row() const
{
    return row_;
}

size_t mind::IndexVector::column() const
{
    return col_;
}

size_t mind::IndexVector::slice() const
{
    return slice_;
}
