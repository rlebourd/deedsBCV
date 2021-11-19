#ifndef Image_hpp
#define Image_hpp

#include <stdio.h>
#include <string>

namespace mind {

    class MRIImage {
    public:
        using NIFTIDimensions = std::tuple<size_t, size_t, size_t, size_t>;
        
        MRIImage(std::string filename);
        
        NIFTIDimensions dimensions() const;
        
        void addScalarToEachVoxel(double scalar);
        
    private:
        std::string filename;
        NIFTIDimensions dim;
    };

}

#endif /* Image_hpp */
