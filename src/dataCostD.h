#ifndef DATA_COST_D_H_INCLUDED
#define DATA_COST_D_H_INCLUDED

void interp3xyz(float* datai,float* data,float* datax,float* datay,int len1,int len2);

void interp3xyzB(float* datai,float* data,float* datax,float* datay,int len1,int len2);

void dataCostCL(unsigned long* data,
                unsigned long* data2,
                float* results,
                int m, int n, int o,
                int len2,
                int step1,
                int hw,
                float quant,
                float alpha,
                int randnum);

void warpImageCL(float* warped,float* im1,float* im1b,float* u1,float* v1,float* w1);

void warpAffineS(short* warped,short* input,float* X,float* u1,float* v1,float* w1);

void warpAffine(float *warped, const float *input, const float *im1b, const float *X, const float *u1, const float *v1, const float *w1);

#endif // DATA_COST_D_H_INCLUDED
