//
//  Image.hpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#ifndef Image_hpp
#define Image_hpp

#include <stdio.h>
#include <string>
#include "Matrix.hpp"

namespace mind {

    /// @brief Represents an immutable 3-dimensional image of black-and-white intensity values
    class Image {
        mind::Matrix intensities_;
        
    public:
        /// @brief Loads the image from a NITFI file with extension nii.gz
        static mind::Image loadedFromFile(std::string filename);
        
        /// @brief Returns a const reference to the image's matrix of intensity values
        const mind::Matrix &intensitiesMatrix() const;
    };

}

#endif /* Image_hpp */
