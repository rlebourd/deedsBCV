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
            for(int i = 0; i < (hw + 1); i++){
                temp2[ind(i, j, k)] = temp1[ind(i+hw, j, k)];
            }
            
            // between the top and bottom borders of temp2
            //
            for(int i = (hw + 1); i < (m - hw); i++){
                temp2[ind(i, j, k)] = temp1[ind(i+hw, j, k)] - temp1[ind(i-hw-1, j, k)];
            }
            
            // on the bottom border of temp2
            //
            for(int i = (m - hw); i < m; i++){
                temp2[ind(i, j, k)] = temp1[ind(m-1, j, k)] - temp1[ind(i-hw-1, j, k)];
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
                temp1[i + j*m + k*m*n] = temp2[i + (j+hw)*m + k*m*n];
            }
            for(int j = (hw + 1); j < (n - hw); j++){
                temp1[i + j*m + k*m*n] = temp2[i + (j+hw)*m + k*m*n] - temp2[i + (j-hw-1)*m + k*m*n];
            }
            for(int j = (n - hw); j < n; j++){
                temp1[i + j*m + k*m*n] = temp2[i + (n-1)*m + k*m*n] - temp2[i + (j-hw-1)*m + k*m*n];
            }
        }
    }
    
    // this integrates along the direction of increasing slice number
    for(int k = 1; k < o; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < m; i++){
                temp1[i + j*m + k*m*n] += temp1[i + j*m + (k-1)*m*n];
            }
        }
    }
    
    //
    for(int j = 0; j < n; j++){
        for(int i = 0; i < m; i++){
            for(int k = 0; k < (hw+1); k++){
                input[i + j*m + k*m*n] = temp1[i+j*m+(k+hw)*m*n];
            }
            for(int k = (hw+1); k < (o-hw); k++){
                input[i + j*m + k*m*n] = temp1[i+j*m+(k+hw)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
            for(int k = (o-hw); k < o; k++){
                input[i + j*m + k*m*n] = temp1[i+j*m+(o-1)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
        }
    }
    
    delete[] temp1; delete[] temp2;
}


void imshift(const float* input,float* output,int dx,int dy,int dz,int m,int n,int o){
    // shifts the image by the given 3D displacement
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                if(i+dy>=0&&i+dy<m&&j+dx>=0&&j+dx<n&&k+dz>=0&&k+dz<o)
                    output[i+j*m+k*m*n]=input[i+dy+(j+dx)*m+(k+dz)*m*n];
                else
                    output[i+j*m+k*m*n]=input[i+j*m+k*m*n];
            }
        }
    }
}

void distances(const float* im1, float* d1, int m, int n, int o, int qs){
    const int sz1 = m*n*o;
    
    float * const w1 = new float[sz1];
    
    // shifts the whole image by a vectorial displacement determined by the input l, which ranges from 0 to 5
    const int dx[6] = {qs,  qs, -qs,   0, qs,  0};
    const int dy[6] = {qs, -qs,   0, -qs,  0, qs};
    const int dz[6] = {0,    0,  qs,  qs, qs, qs};
    
#pragma omp parallel for
    for(int l=0; l < 6; l++){
        // shifts the image by dx[l], dy[l], dz[l]
        imshift(im1, w1, dx[l], dy[l], dz[l], m, n, o);
        
        // calculates the squared distance between each voxel in the original image and the shifted image
        for(int i = 0; i < sz1; i++){
            const auto dw = w1[i] - im1[i];
            w1[i] = dw*dw;
        }
        
        // blurs the squared-distances matrix using qs as the blur radius
        boxfilter(w1, qs, m, n, o);
        for(int i = 0; i < sz1; i++){
            // stores the whole squared-distances matrix in level-l of the 4-D matrix d1
            d1[i + l*sz1] = w1[i];
            
            // stores 6 floats per voxel => d1[i + j*m + k*m*n + l*m*n*o] =
            // component l of the vector for voxel (i, j, k)
            // so we can see here that d1 is 4-dimensional and it is
            // accessed like d1(i, j, k, l).
        }
    }
    
    delete[] w1;
}

//__builtin_popcountll(left[i]^right[i]); absolute hamming distances
void descriptor(uint64_t* mindq, float* im1, int m, int n, int o, int qs){
    //MIND with self-similarity context
    const auto ind = [&](int i, int j, int k){
        return i + j*m + k*m*n;
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
    const int numberOfVoxels = m*n*o;
    const int numberOfDistancesPerVoxel = 6; // (aka the number of vectorial displacements r in the search region R)
    float* d1 = new float[numberOfVoxels*numberOfDistancesPerVoxel]; // 6 floats per voxel
    distances(im1, d1, m, n, o, qs);
    
    //
    const int sx[12] = {-qs,   0, -qs,  0,   0, qs,   0,  0,   0, -qs,   0,   0};
    const int sy[12] = {  0, -qs,   0, qs,   0,  0,   0, qs,   0,   0,   0, -qs};
    const int sz[12] = {  0,   0,   0,  0, -qs,  0, -qs,  0, -qs,   0, -qs,   0};
    
    const int index[12] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
    const int len2 = 12;
    image_d = 12;
    
    //quantisation table
    const int val = 6;
    
#pragma omp parallel for
    float mindThreshold[val-1];
    for(int i = 0; i < val-1; i++){
        mindThreshold[i] = -log((i+1.5f)/val);
    }
    
    for(int k = 0; k < o; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < m; i++){
                float mind1[12];

                // initialize the mind descriptors from the blurred distances
                for(int l = 0; l < len2; l++){
                    if((i + sy[l]) >= 0 &&
                       (i + sy[l]) < m &&
                       (j + sx[l]) >= 0 &&
                       (j + sx[l]) < n &&
                       (k + sz[l]) >= 0 &&
                       (k + sz[l]) < o)
                    {
                        mind1[l] = d1[i + sy[l] + (j+sx[l])*m + (k+sz[l])*m*n + index[l]*numberOfVoxels];
                    }
                    else{
                        mind1[l] = d1[ind(i,j,k) + index[l]*numberOfVoxels];
                    }
                }
                
                // calculate the Z-score of the mind descriptor
                // (value - mean) / std-dev
                const float minimumBlurredSquaredDistance = *min_element(mind1, mind1 + len2);
                float totalNoise = 0.0f;
                for(int l = 0; l < len2; l++){
                    mind1[l] -= minimumBlurredSquaredDistance;
                    totalNoise += mind1[l];
                }
                const float averageNoise = max(totalNoise/(float)len2, 1e-6f);
                for(int l = 0; l < len2; l++){
                    mind1[l] /= averageNoise;
                }
                
                //
                unsigned long long accum=0;
                unsigned long long tabled1=1;
                
                for(int l = 0; l < len2; l++){
                    const unsigned long long power = 32;
                    const unsigned int tablei[6] = {0, 1, 3, 7, 15, 31};

                    int mind1val = 0;
                    for(int c = 0; c < val-1; c++){
                        mind1val += (mindThreshold[c] > mind1[l]) ? 1 : 0;
                    }
                    accum += tablei[mind1val] * tabled1;
                    tabled1 *= power;
                }
                mindq[ind(i,j,k)] = accum;
            }
        }
    }
    
    delete[] d1;
}


