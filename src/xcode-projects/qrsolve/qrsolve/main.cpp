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

float weightedAverageOfPixelAt(lin_alg::Vector4D xyz, mind::MRIImage input){
    const float x1=xyz.at(0);
    const float y1=xyz.at(1);
    const float z1=xyz.at(2);
    
    const int x=floor(x1);
    const int y=floor(y1);
    const int z=floor(z1);

    const float dx=x1-x; // in range [0, 1] and measures distance of index from nearest (lower bound) whole number
    const float dy=y1-y; // in range [0, 1] so if dy ~ 0 then (1-dy) ~ 1, and y is closer to floor(y) than to floor(y)+1
    const float dz=z1-z; // in range [0, 1] so if dz ~ 1 then (1-dz) ~ 0, and z is closer to floor(z)+1 than to floor(z)

    return lin_alg::Matrix2x2{{
        input.at(x+0, y+0, z+0), input.at(x+0, y+1, z+0),
        input.at(x+1, y+0, z+0), input.at(x+1, y+1, z+0)
    }}.dottedWith(lin_alg::Matrix2x2{{
        (1.0-dx)*(1.0-dy)*(1.0-dz), (1.0-dx)*dy*(1.0-dz),
        dx*(1.0-dy)*(1.0-dz), dx*dy*(1.0-dz)
    }})
    
    +
    
    lin_alg::Matrix2x2{{
        input.at(x+0, y+0, z+1), input.at(x+0, y+1, z+1),
        input.at(x+1, y+0, z+1), input.at(x+1, y+1, z+1)
    }}.dottedWith(lin_alg::Matrix2x2{{
        (1.0-dx)*(1.0-dy)*dz, (1.0-dx)*dy*dz,
        dx*(1.0-dy)*dz, dx*dy*dz
    }});
}

mind::MRIImage warpAffine(mind::MRIImage deformedImage, lin_alg::Matrix4x4 affineMatrix){
    auto originalImage = deformedImage;
    
    const int o = static_cast<int>(std::get<2>(originalImage.dimensions()));
    const int n = static_cast<int>(std::get<1>(originalImage.dimensions()));
    const int m = static_cast<int>(std::get<0>(originalImage.dimensions()));

    for (int k=0; k<o; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<m; i++){
                // visit each voxel in the undeformed image (i.e. indexed by i, j, k)
                
                // multiply the column vector of input indices by the affine
                // transformation matrix to calculate the corresponding indices
                // in the deformed image
                const auto originalIndices = lin_alg::Vector4D{{(double)i, (double)j, (double)k, 1}};
                const auto deformedIndices = affineMatrix.times(originalIndices);
                
                // would be nice to replace the if-else block below with:
                // originalImage.at(originalIndices) = deformedImage.at(deformedIndices)
                //
                // i don't like this because it isolates the voxel-averaging algorithm
                // within the image class, which seems inappropriate as it's not extensible
                
                // or something like
                // originalImage.at(originalIndices) = deformedImage.reduceForeachVoxelNear(0, deformedIndices, [](const VoxelProxy &vox){
                //    return vox.value * vox.dx * vox.dy * vox.dz;
                // });
                //
                // i like this because it's abstract, intuitive, and extensible
                
                if (deformedImage.pixelsAreOnOrOutsideBoundary(deformedIndices)){
                    // ignore indices on the boundary of the image
                    originalImage.at(i, j, k) = 0.0;
                } else {
                    // deform the undeformed image
                    
                    // this is just computing the weighted average of the pixels near the real-valued index
                    originalImage.at(i, j, k) = weightedAverageOfPixelAt(deformedIndices, deformedImage);
                }
            }
        }
    }
    return originalImage;
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
