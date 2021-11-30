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

int main(int argc, const char * argv[]) {
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
            const auto F1 = integratedDifferences.at(idx.subtracting(patchVector));
            const auto F2 = integratedDifferences.at(idx.adding(patchVector));
            return (F2 - F1);
        };
        distances[dr] = integratedDifferences.mappingElements(convolution);
    }
    
    return 0;
}
