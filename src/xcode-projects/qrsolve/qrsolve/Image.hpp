#ifndef Image_hpp
#define Image_hpp

#include <stdio.h>
#include <string>
#include <vector>

namespace mind {

    class MRIImage {
    public:
        using NIFTIDimensions = std::tuple<size_t, size_t, size_t, size_t>;
        
        MRIImage(std::string filename);
        
        NIFTIDimensions dimensions() const;
        
        void addScalarToVoxels(double scalar);
        
        double at(size_t i, size_t j, size_t k) const;

        double &at(size_t i, size_t j, size_t k);
        
    private:
        std::string filename;
        NIFTIDimensions dimensionsOfImage;
        std::vector<double> imageData;
    };

}

#endif /* Image_hpp */
