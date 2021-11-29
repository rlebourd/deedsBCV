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
        
    };
    
    const size_t patchLength = 5u;
    const auto patchVector = mind::IndexVector{patchLength, patchLength, patchLength};
    const auto originalImage = mind::Image::loadedFromNitfiFile("example.nii.gz");
        
    const auto &originalIntensities = originalImage.intensitiesMatrix();
    for (auto dr : searchRegion){
        const mind::Matrix integratedDifferences = originalIntensities
                                                       .subtracting(originalIntensities.translatedBy(dr))
                                                       .integratedAlongXYZ();
        distances[dr] = integratedDifferences.mappingElements([](mind::IndexVector &idx, double element){
            // TODO: calculate the convolution
            
        });
    }
    
    return 0;
}
