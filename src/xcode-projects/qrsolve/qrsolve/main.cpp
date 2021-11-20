//
//  main.cpp
//  qrsolve
//
//  Created by Robert LeBourdais on 11/18/21.
//

#include <iostream>
#include "Matrix.hpp"
#include "Image.hpp"
#include <cassert>
#include <algorithm>
#include <math.h>

lin_alg::Matrix4x4 qrsolve(lin_alg::Matrix4x4 &Ainv, const lin_alg::Matrix4x4 &A){
    return Ainv;
}

mind::MRIImage warpAffine(mind::MRIImage input, lin_alg::Matrix4x4 affineMatrix){
    auto output = input;
    
    const int o = static_cast<int>(std::get<2>(output.dimensions()));
    const int n = static_cast<int>(std::get<1>(output.dimensions()));
    const int m = static_cast<int>(std::get<0>(output.dimensions()));

    for(int k=0; k<o; k++){
        for(int j=0; j<n; j++){
            for(int i=0; i<m; i++){
                // visit each voxel in the undeformed image (i.e. indexed by i, j, k)
                
                // multiply the column vector of input indices by the affine
                // transformation matrix to calculate the corresponding indices
                // in the deformed image
                auto originalIndices = lin_alg::Vector4D{{static_cast<double>(i),
                                                          static_cast<double>(j),
                                                          static_cast<double>(k),
                                                          1}};
                auto deformedIndices = affineMatrix.times(originalIndices);

                const float x1=deformedIndices.at(0);
                const float y1=deformedIndices.at(1);
                const float z1=deformedIndices.at(2);
                
                const int x=floor(x1);
                const int y=floor(y1);
                const int z=floor(z1);
                
                if(y<0 || y>=m-1 || x<0 || x>=n-1 || z<0 || z>=o-1){
                    // ignore indices on the boundary of the image
                    output.at(i, j, k) = 0.0;
                }
                else{
                    // deform the undeformed image
                    const float dx=x1-x; // in range [0, 1] and measures distance of index from nearest (lower bound) whole number
                    const float dy=y1-y; // in range [0, 1] so if dy ~ 0 then (1-dy) ~ 1, and y is closer to floor(y) than to floor(y)+1
                    const float dz=z1-z; // in range [0, 1] so if dz ~ 1 then (1-dz) ~ 0, and z is closer to floor(z)+1 than to floor(z)
                    
                    // this is just computing the weighted average of the pixels near the real-valued index
                    output.at(i, j, k) =
                        input.at(x+0, y+0, z+0) * (1.0-dx)*(1.0-dy)*(1.0-dz) +
                        input.at(x+0, y+1, z+0) * (1.0-dx)*dy*(1.0-dz) +
                        input.at(x+1, y+0, z+0) * dx*(1.0-dy)*(1.0-dz) +
                        input.at(x+0, y+0, z+1) * (1.0-dx)*(1.0-dy)*dz +
                        input.at(x+1, y+1, z+0) * dx*dy*(1.0-dz) +
                        input.at(x+0, y+1, z+1) * (1.0-dx)*dy*dz +
                        input.at(x+1, y+0, z+1) * dx*(1.0-dy)*dz +
                        input.at(x+1, y+1, z+1) * dx*dy*dz;
                }
            }
        }
    }
    return output;
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
    
#if 0
    // the initialization of the linear BCV algorithm goes something like
    auto movingImage = mind::MRIImage{"moving.nii.gz"};
    auto fixedImage  = mind::MRIImage{"fixed.nii.gz"};
    if (movingImage.dimensions() != fixedImage.dimensions()){
        return -1;
    }
    
    movingImage.addScalarToEachVoxel(1024);
    fixedImage.addScalarToEachVoxel(1024);

    const auto X = lin_alg::Matrix4x4::identityMatrix();
    const auto Xprev = X;
#endif // 0
    
    
    return 0;
}
