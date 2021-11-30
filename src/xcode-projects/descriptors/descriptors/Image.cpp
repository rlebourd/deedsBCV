//
//  Image.cpp
//  descriptors
//
//  Created by Robert LeBourdais on 11/29/21.
//

#include "Image.hpp"

mind::Image::Image(mind::Matrix intensities)
: intensities_(intensities)
{
    
}

mind::Image mind::Image::loadedFromFile(std::string filename)
{
    const size_t nrows = 0;
    const size_t ncols = 0;
    const size_t nslices = 0;
    const size_t nelements = nrows*ncols*nslices;
    std::vector<double> elems = std::vector<double>(0, nelements);
    
    // TODO: open the file
    // TODO: set the nrows, ncols, and nslices from the file data
    // TODO: populate the elems vector from the file data
    
    mind::Matrix intensities = mind::Matrix::fromData(nrows, ncols, nslices, elems);
    mind::Image image{intensities};
    return image;
}

const mind::Matrix &mind::Image::intensitiesMatrix() const
{
    return intensities_;
}
