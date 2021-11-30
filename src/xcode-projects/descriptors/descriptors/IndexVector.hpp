//
//  Vector.hpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#ifndef Vector_hpp
#define Vector_hpp

#include <stdio.h>

namespace mind {

    class IndexVector {
        size_t row_;
        size_t col_;
        size_t slice_;
        
    public:
        IndexVector(size_t row, size_t col, size_t slice);
        
        IndexVector subtracting(const IndexVector &rhs) const;
        IndexVector adding(const IndexVector &rhs) const;
        size_t row() const;
        size_t column() const;
        size_t slice() const;
    };

}

#endif /* Vector_hpp */
