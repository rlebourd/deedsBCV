#ifndef MIND_SSC_BOX_H_INCLUDED
#define MIND_SSC_BOX_H_INCLUDED

#include <stdint.h>

void boxfilter(float* input, int hw, int m, int n, int o);

void imshift(const float* input,float* output,int dx,int dy,int dz,int m,int n,int o);

void distances(const float* im1,float* d1,int m,int n,int o,int qs);

void descriptor(uint64_t* mindq,float* im1,int m,int n,int o,int qs);

#endif // MIND_SSC_BOX_H_INCLUDED
