//
//  main.cpp
//  qrsolve
//
//  Created by Robert LeBourdais on 11/18/21.
//

#include <iostream>
#include "Matrix.hpp"
#include <cassert>

lin_alg::Matrix4x4 qrsolve(lin_alg::Matrix4x4 &Ainv, const lin_alg::Matrix4x4 &A){
    
    
    return Ainv;
}

int main(int argc, const char * argv[]) {
    const auto I  = lin_alg::Matrix3x3::identityMatrix();
    const auto I2 = I.times(2);
    const auto I4 = I.times(4);
    const auto A = I2.times(I2);
    assert(A == I4);
    
    const auto x = lin_alg::Vector3D{{1, 2, 3}};
    const auto b = A.times(x);
    assert((b == lin_alg::Vector3D{{4, 8, 12}}));
    
    return 0;
}
