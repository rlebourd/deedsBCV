#include "MINDSSCbox.h"
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>

extern bool RIGID;
extern int RAND_SAMPLES;
extern int image_m;
extern int image_n;
extern int image_o;
extern int image_d;
extern float SSD0;
extern float SSD1;
extern float SSD2;
extern float distfx_global;
extern float beta;
extern int qc;

using namespace std;

void boxfilter(float* input, int hw, int m, int n, int o){
    const auto ind = [&](int i, int j, int k){
        return i + j*m + k*m*n;
    };

    const int sz = m*n*o;
    float * const temp1 = new float[sz];
    float * const temp2 = new float[sz];

    // this copies the input image to the temporary buffer
    for(int i=0; i < sz; i++){
        temp1[i] = input[i];
    }
    
    // this integrates the image along the direction of increasing row number
    for(int k = 0; k < o; k++){
        for(int j = 0; j < n;j++){
            for(int i = 1; i < m ;i++){
                temp1[ind(i, j, k)] += temp1[ind(i-1, j, k)];
            }
        }
    }
    
    //
    for(int k = 0; k < o; k++){
        for(int j = 0; j < n; j++){
            // on the top border of temp2
            //
            // This copies the integral of the image from the rows in the range [hw, 2*hw]
            // into the rows in the range [0, hw]. Why is this necessary?
            //
            for(int i = 0; i < (hw + 1); i++){
                temp2[ind(i, j, k)] = temp1[ind(i+hw, j, k)]; // note: row (i-hw-1) < (hw+1)-hw-1 = 0 => that row does not exist yet
            }
            
            // between the top and bottom borders of temp2
            //
            for(int i = (hw + 1); i < (m - hw); i++){
                temp2[ind(i, j, k)] = temp1[ind(i+hw, j, k)] - temp1[ind(i-hw-1, j, k)];
            }
            
            // on the bottom border of temp2
            //
            for(int i = (m - hw); i < m; i++){
                temp2[ind(i, j, k)] = temp1[ind(m-1, j, k)] - temp1[ind(i-hw-1, j, k)]; // note: last row minus (hw+1) before current row
                // this reflects an attempt to "sort of integrate"
                // the intensity to the bottom row
            }
        }
    }
    
    // this integrates temp2 along the direction of increasing column number
    for(int k = 0; k < o; k++){
        for(int j = 1; j < n; j++){
            for(int i = 0; i < m; i++){
                temp2[ind(i, j, k)] += temp2[ind(i, j-1, k)];
            }
        }
    }
    
    //
    for(int k = 0; k < o; k++){
        for(int i=0; i < m; i++){
            for(int j = 0; j < (hw + 1); j++){
                temp1[ind(i, j, k)] = temp2[ind(i, j+hw, k)];
            }
            for(int j = (hw + 1); j < (n - hw); j++){
                temp1[ind(i, j, k)] = temp2[ind(i, j+hw, k)] - temp2[ind(i, j-hw-1, k)];
            }
            for(int j = (n - hw); j < n; j++){
                temp1[ind(i, j, k)] = temp2[ind(i, n-1, k)] - temp2[ind(i, j-hw-1, k)];
            }
        }
    }
    
    // this integrates along the direction of increasing slice number
    for(int k = 1; k < o; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < m; i++){
                temp1[ind(i, j, k)] += temp1[ind(i, j, k-1)];
            }
        }
    }
    
    //
    for(int j = 0; j < n; j++){
        for(int i = 0; i < m; i++){
            for(int k = 0; k < (hw+1); k++){
                input[ind(i, j, k)] = temp1[ind(i, j, k+hw)];
            }
            for(int k = (hw+1); k < (o-hw); k++){
                input[ind(i, j, k)] = temp1[ind(i, j, k+hw)]-temp1[ind(i, j, k-hw-1)];
            }
            for(int k = (o-hw); k < o; k++){
                input[ind(i, j, k)] = temp1[ind(i, j, o-1)]-temp1[ind(i, j, k-hw-1)];
            }
        }
    }
    
    delete[] temp1; delete[] temp2;
}


void shiftImage(const float* input,float* output,int dx,int dy,int dz,int m,int n,int o){
    // shifts the image by the given 3D displacement
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                if(i+dy>=0&&i+dy<m&&j+dx>=0&&j+dx<n&&k+dz>=0&&k+dz<o)
                    output[i+j*m+k*m*n] = input[i+dy+(j+dx)*m+(k+dz)*m*n];
                else
                    output[i+j*m+k*m*n] = input[i+j*m+k*m*n];
            }
        }
    }
}

void distances(const float* originalImage, float** d1, int m, int n, int o, int qs){
    const int numberOfVoxelsInImage = m*n*o;
    const int numberOfPositionsInSearchRegion = 6; // (aka the number of vectorial displacements r in the search region R)
    *d1 = new float[numberOfVoxelsInImage*numberOfPositionsInSearchRegion]; // 6 floats per voxel

    float * const squaredDistancesForSearchPosition = new float[numberOfVoxelsInImage];
    
    // shifts the whole image by a vectorial displacement
    // determined by the input l, which ranges from 0 to 5
    //
    // this traverses the "search region," which consists of
    // 6 points in the neighborhood of the central voxel.
    //
    const int searchPosition_dx[numberOfPositionsInSearchRegion] = {qs,  qs, -qs,   0, qs,  0};
    const int searchPosition_dy[numberOfPositionsInSearchRegion] = {qs, -qs,   0, -qs,  0, qs};
    const int searchPosition_dz[numberOfPositionsInSearchRegion] = {0,    0,  qs,  qs, qs, qs};
    
#pragma omp parallel for
    for(int searchPositionIndex = 0; searchPositionIndex < numberOfPositionsInSearchRegion; searchPositionIndex++){
        // shifts the image by dx, dy, dz
        const auto dx = searchPosition_dx[searchPositionIndex];
        const auto dy = searchPosition_dy[searchPositionIndex];
        const auto dz = searchPosition_dz[searchPositionIndex];
        shiftImage(originalImage,
                   squaredDistancesForSearchPosition,
                   dx, dy, dz,
                   m, n, o);
        
        // calculates the squared distance between each voxel in the original image and the shifted image
        for(int i = 0; i < numberOfVoxelsInImage; i++){
            const auto dw = squaredDistancesForSearchPosition[i] - originalImage[i];
            squaredDistancesForSearchPosition[i] = dw*dw;
        }
        
        // blurs the squared-distances matrix using qs as the blur radius
        boxfilter(squaredDistancesForSearchPosition, qs, m, n, o);
        for(int i = 0; i < numberOfVoxelsInImage; i++){
            // stores the whole squared-distances matrix at depth "searchPositionIndex" of the 4-D matrix d1
            (*d1)[i + searchPositionIndex*numberOfVoxelsInImage] = squaredDistancesForSearchPosition[i];
            
            // stores 6 floats per voxel => d1[i + j*m + k*m*n + l*m*n*o] =
            // component l of the vector for voxel (i, j, k)
            // so we can see here that d1 is 4-dimensional and it is
            // accessed like d1(i, j, k, l).
        }
    }
    
    delete[] squaredDistancesForSearchPosition;
}

//__builtin_popcountll(left[i]^right[i]); absolute hamming distances
void descriptor(uint64_t* mindq, float* im1, int m, int n, int o, int qs){
    //MIND with self-similarity context
    const auto ind3 = [&](int i, int j, int k){
        return i + j*m + k*m*n;
    };
    const auto ind4 = [&](int i, int j, int k, int l){
        return i + j*m + k*m*n + l*m*n*o;
    };

    //============== DISTANCES USING BOXFILTER ===================
    
    // Calculate blurred distances for 6 different vectorial shifts of the image
    //
    // This calculation yields the distances function D(I, x, x+r) that returns the
    // Euclidean distance between the image patches of size p centered at positions
    // x and x+r in R3.
    //
    // The size p of the patch is determined by the parameter qs, which is calculated
    // as floor(quantisation/2 + 1) and which is called the "mind step" elsewhere in
    // the program.
    //
    // When we say that that "patch size is p," what we mean is that the side length
    // of the patch is 2p+1. That is, there are p voxels directly to the left and p
    // voxels directly to the right of the voxel at the center of the patch, making
    // the side length equal to 2p+1. Similar statements can be made for any of the
    // three coordinate axes.
    //
    float *d1 = nullptr;
    distances(im1, &d1, m, n, o, qs);
    
    // A number of shifts
    const int sx[12] = {-qs,   0, -qs,  0,   0, qs,   0,  0,   0, -qs,   0,   0};
    const int sy[12] = {  0, -qs,   0, qs,   0,  0,   0, qs,   0,   0,   0, -qs};
    const int sz[12] = {  0,   0,   0,  0, -qs,  0, -qs,  0, -qs,   0, -qs,   0};
    
#if 0
    // compare these shifts to:
    
    const int searchPosition_dx[numberOfPositionsInSearchRegion] = {qs,  qs, -qs,   0, qs,  0};
    const int searchPosition_dy[numberOfPositionsInSearchRegion] = {qs, -qs,   0, -qs,  0, qs};
    const int searchPosition_dz[numberOfPositionsInSearchRegion] = {0,    0,  qs,  qs, qs, qs};

    dx[0] = -(sx[0]+sx[1])
    dy[0] = -(sy[0]+sy[1])
    dz[0] = -(sz[0]+sz[1])

    dx[1] = -(sx[2]+sx[3])
    dy[1] = -(sy[2]+sy[3])
    dz[1] = -(sz[2]+sz[3])

    dx[2] = -(sx[4]+sx[5])
    dy[2] = -(sy[4]+sy[5])
    dz[2] = -(sz[4]+sz[5])

    // more generally:
    
    dx[k] = -(sx[2*k] + sx[2*k+1])
    dy[k] = -(sy[2*k] + sy[2*k+1])
    dz[k] = -(sz[2*k] + sz[2*k+1])
    
    // why this happens to be the case is not yet clear.

#endif // 0
    
    const int shiftIndexToSearchRegionIndex[12] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
    const int numberOfNeighbors = 12;
    image_d = 12;
    
    // Quantisation table
    //
    // The original paper on MIND states that,
    //
    //   "In order to speed-up these computations
    //    the descriptor can be quantised to only
    //    4 bit, without significant loss of accuracy.
    //
    //    For |R| = 6 all possible distances between
    //    descriptors can be pre-computed and stored
    //    in a lookup-table."
    //
    // A 4 bit number can represent the numbers 0 to 15
    // in decimal, 0 to 0xF in hexadecimal, and 0 to 1111
    // in binary. So it can represent 16 different numbers.
    //
    // If |R| = 6, then the MIND descriptor has six 4-bit
    // components, so each component can take on one of
    // 16 different values. There are 16^6 different
    // descriptors having six 4-bit components. This is
    // 16.7 million possible descriptors, with each
    // descriptor occupying 24 bits (3 bytes) of memory.
    // Therefore, a table of all such descriptors would
    // take up about 50 million bytes of memory, or 50 MB,
    // which is a modest amount of memory even on a
    // smartphone.
    //
    // Since there are 16.7 million different possible
    // descriptors, we cannot store a lookup table
    // that contains the distance between any two of
    // these 16.7 million descriptors, because that would
    // require about 280 trillion pairs of descriptors,
    // which would far exceed the memory capacity of any
    // modern computer.
    //
    // Perhaps when he says that the descriptor can be
    // quantised to 4 bit, he means the entire descriptor
    // and not just a single component. In this interpretation,
    // it would be possible to store all possible distances
    // between descriptors, because there would only be 256
    // different pairs of descriptors. But that doesn't
    // seem to make sense as an interpretation because the
    // datatype for a descriptor is uint64_t.
    //
    // Perhaps what he means is that the distance is the
    // distance between individual MIND *components* not
    // between MIND *vectors*. This would reconcile all
    // the available facts, in that the whole descriptor
    // would require 24 bits of RAM, meaning that the
    // uint64_t datatype would make some more sense as a
    // choice, and that the distance between any pair of
    // corresponding MIND *components* could be pre-
    // computed and stored in a lookup table, specifically
    // because there are 256 such possible pairs.
    //
    // With that said, I don't see a variable that stands
    // out as a lookup table for precomputed distances,
    // despite the fact that the original code was labeled
    // as "Quantisation table."
    //
    const int numberOfMindThresholds = 6;
    float mindThreshold[numberOfMindThresholds];
    for(int i = 0; i < numberOfMindThresholds; i++){
        mindThreshold[i] = -log((i+1.5f)/(numberOfMindThresholds+1));
    }
    
#pragma omp parallel for
    for(int k = 0; k < o; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < m; i++){
                float distancesForCurrentVoxel[numberOfNeighbors];

                // initialize the mind descriptors from the blurred distances
                for(int shiftIndex = 0; shiftIndex < numberOfNeighbors; shiftIndex++){
                    const int searchRegionIndex = shiftIndexToSearchRegionIndex[shiftIndex];
                    
                    if((i + sy[shiftIndex]) >= 0 &&
                       (i + sy[shiftIndex]) <  m &&
                       (j + sx[shiftIndex]) >= 0 &&
                       (j + sx[shiftIndex]) <  n &&
                       (k + sz[shiftIndex]) >= 0 &&
                       (k + sz[shiftIndex]) <  o)
                    {
                        // the shifted pixel is in bounds
                        //
                        // this is the typical execution path,
                        // except for voxels on the boundary
                        // of the 3D image
                        //
                        //
                        distancesForCurrentVoxel[shiftIndex] = d1[ind4(i+sy[shiftIndex],
                                                                       j+sx[shiftIndex],
                                                                       k+sz[shiftIndex],
                                                                       searchRegionIndex)];
                    }
                    else{
                        // the shifted pixel is out of bounds
                        //
                        // this case only executes for voxels
                        // on the boundary of the 3D image
                        //
                        // so this is an edge case.
                        //
                        distancesForCurrentVoxel[shiftIndex] = d1[ind4(i,
                                                                       j,
                                                                       k,
                                                                       searchRegionIndex)];
                    }
                }
                
                // calculate the Z-score of the mind descriptor
                // (value - mean) / std-dev
                const float minimumDistanceFromVoxelToNeighbor = *min_element(distancesForCurrentVoxel, distancesForCurrentVoxel + numberOfNeighbors);
                float totalNoise = 0.0f;
                for(int shiftIndex = 0; shiftIndex < numberOfNeighbors; shiftIndex++){
                    distancesForCurrentVoxel[shiftIndex] -= minimumDistanceFromVoxelToNeighbor;
                    totalNoise += distancesForCurrentVoxel[shiftIndex];
                }
                const float averageNoise = max(totalNoise/(float)numberOfNeighbors, 1e-6f);
                for(int neighborIndex = 0; neighborIndex < numberOfNeighbors; neighborIndex++){
                    distancesForCurrentVoxel[neighborIndex] /= averageNoise;
                }
                
                //
                // A reading of the software below seems to reveal that the MIND
                // descriptor has *twelve*, and NOT six, 4-bit components.
                //
                // This can be seen by observing that, since numberOfShifts = 12,
                // tabled1, which is initially equal to 1, is multiplied by 32
                // twelve times. But (2^5)^12 = 2^60, which means that the
                // the MIND descriptor requires 60 bits of memory in this
                // representation. This explains why it is stored as a 64-bit
                // value.
                unsigned long long accum=0;
                unsigned long long tabled1=1;
                
                for(int shiftIndex = 0; shiftIndex < numberOfNeighbors; shiftIndex++){
                    const unsigned long long power = 32;
                    const unsigned int tablei[6] = {0, 1, 3, 7, 15, 31};

                    int mind1val = 0;
                    for(int threshIndex = 0; threshIndex < numberOfMindThresholds; threshIndex++){
                        mind1val += (mindThreshold[threshIndex] > distancesForCurrentVoxel[shiftIndex]) ? 1 : 0;
                    }
                    accum += tablei[mind1val] * tabled1;
                    tabled1 *= power;
                }
                mindq[ind3(i,j,k)] = accum;
            }
        }
    }
    
    delete[] d1;
}


