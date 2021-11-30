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
    
    const auto searchRegion = std::vector<mind::IndexVector>{
        mind::IndexVector{1, 0, 0},
        mind::IndexVector{0, 1, 0},
        mind::IndexVector{0, 0, 1}
    };
    
    const size_t patchLength = 5u;
    const auto patchVector = mind::IndexVector{patchLength, patchLength, patchLength};
    const auto originalImage = mind::Image::loadedFromFile("example.nii.gz");
    
    const auto &originalIntensities = originalImage.intensitiesMatrix();
    for (auto dr : searchRegion){
        const mind::Matrix integratedDifferences = originalIntensities
                                                       .subtracting(originalIntensities.translatedBy(dr))
                                                       .integratedAlongXYZ();
        const auto convolution = [&](mind::IndexVector &idx, double element){
            // TODO: address underflow and overflow
            //
            // adding the patch vector could put the index
            // beyond the bounds of the matrix subtracting
            // the patch vector could put the index below
            // (0, 0, 0), also outside the bounds of the
            // matrix.
            const auto F1 = integratedDifferences.at(idx.subtracting(patchVector));
            const auto F2 = integratedDifferences.at(idx.adding(patchVector));
            return (F2 - F1);
        };
        distances[dr] = integratedDifferences.mappingOverElements(convolution);
    }
    
    return distances;
}

int main(int argc, const char * argv[]) {
    auto distances = createDistances();
    
    return 0;
}
