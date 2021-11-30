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
    // create the search region
    const size_t qs = 0u;
    const auto searchRegion = std::vector<mind::IndexVector>{
        mind::IndexVector{ qs,  qs,  0},
        mind::IndexVector{ qs, -qs,  0},
        mind::IndexVector{-qs,   0, qs},
        mind::IndexVector{  0, -qs, qs},
        mind::IndexVector{ qs,   0, qs},
        mind::IndexVector{  0,  qs, qs},
    };
    
    // create the patch vector
    const size_t patchLength = 5u;
    const auto patchVector = mind::IndexVector{patchLength,
                                               patchLength,
                                               patchLength};
    
    // create the image
    const auto originalImage = mind::Image::loadedFromFile("example.nii.gz");
    const auto &originalIntensities = originalImage.intensitiesMatrix();
    
    // populate the distances lookup table
    std::map<mind::IndexVector, mind::Matrix> distances;
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
            //
            // Identify how this error can be handled. One option
            // is to handle the overflow or underflow by defaulting
            // to a reasonable value, but arguably that hides a non-
            // obvious implementation detail from explicit identification.
            //
            // It seems better to me to have the overflow or underflow
            // reported by an exception. If the overflow or underflow
            // occurs, then we default to a reasonable value explicitly.
            //
            // Another alternative is to design the API provided by the
            // matrix class so that it is not possible to produce an
            // exception, but where underflow and overflow are treated
            // explicitly. The API would be something like:
            //
            // idx.subtractingWithDefaultValueOnUnderflow(patchVector, idx);
            // idx.addingWithDefaultValueOnOverflow(patchVector, idx);
            //
            // This models the cases of underflow and overflow explicitly,
            // drawing attention to them as important algorithmic details
            // while providing an expressive API.
            //
            auto F1 = integratedDifferences.at(idx);
            try {
                F1 = integratedDifferences.at(idx.subtracting(patchVector));
            } catch (...){
                // do nothing in case of underflow
            }

            auto F2 = integratedDifferences.at(idx);
            try {
                F2 = integratedDifferences.at(idx.adding(patchVector));
            } catch (...){
                // do nothing in case of overflow
            }
            
            //const auto F1 = integratedDifferences.at(idx.subtracting(patchVector));
            //const auto F2 = integratedDifferences.at(idx.adding(patchVector));
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
