#ifndef TRANSFORMATIONS_H_INCLUDED
#define TRANSFORMATIONS_H_INCLUDED

/* several functions to interpolate and symmetrise deformations
 calculates Jacobian and harmonic Energy */


void interp3(float* interp,float* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag);

void filter1(float* imagein,float* imageout,int m,int n,int o,float* filter,int length,int dim);

void volfilter(float* imagein,int m,int n,int o,int length,float sigma);

float jacobian(float* u1,float* v1,float* w1,int m,int n,int o,int factor);

void consistentMappingCL(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor);

void upsampleDeformationsCL(float* u1,float* v1,float* w1,float* u0,float* v0,float* w0,int m,int n,int o,int m2,int n2,int o2);

#endif // TRANSFORMATIONS_H_INCLUDED
