//
//  main.cpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#include <iostream>
#include "IndexVector.hpp"
#include "Matrix.hpp"
#include "Image.hpp"
#include <map>

static std::map<mind::IndexVector, mind::Matrix> createDistances(){
    std::map<mind::IndexVector, mind::Matrix> distances;
    
    const size_t qs = 0u;
    const auto searchRegion = std::vector<mind::IndexVector>{
        mind::IndexVector{ qs,  qs,  0},
        mind::IndexVector{ qs, -qs,  0},
        mind::IndexVector{-qs,   0, qs},
        mind::IndexVector{  0, -qs, qs},
        mind::IndexVector{ qs,   0, qs},
        mind::IndexVector{  0,  qs, qs},
    };
    
    const size_t patchLength = 5u;
    const auto patchVector = mind::IndexVector{patchLength, patchLength, patchLength};
    const auto originalImage = mind::Image::loadedFromFile("example.nii.gz");
    
    const auto &originalIntensities = originalImage.intensitiesMatrix();
    for (auto dr : searchRegion){
        const mind::Matrix integratedDifferences = originalIntensities
                                                       .subtracting(originalIntensities.translatedBy(dr))
                                                       .integratedAlongXYZ();
        const auto convolution = [&](const mind::IndexVector &idx, double element){
            // TODO: address underflow and overflow
            //
            // Adding the patch vector could put the index
            // beyond the bounds of the matrix, and subtracting
            // the patch vector could put the index below
            // the bounds of the matrix (i.e with negative
            // components)
            const auto F1 = integratedDifferences.at(idx.subtracting(patchVector));
            const auto F2 = integratedDifferences.at(idx.adding(patchVector));
            return (F2 - F1);
        };
        distances[dr] = integratedDifferences.mappingOverElements(convolution);
    }
    
    return distances;
}

int main(int argc, const char * argv[]) {
    // create the SSC-based descriptors
    auto distances = createDistances();
    
    // normalize the distances
    
    
    // quantize the components of the descriptor
    
    
    // produce the descriptor
    
    
    return 0;
}
