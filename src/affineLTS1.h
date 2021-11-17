#ifndef AFFINE_LT1_H_INCLUDED
#define AFFINE_LT1_H_INCLUDED

/// @brief Multiplies a pair of 4x4 matrices
/// @param [in] A is a 4x4 matrix
/// @param [in] B is a 4x4 matrix
/// @param [in] C is a 4x4 matrix
///
/// @note These matrices are represented as arrows of elements in row order, i.e. {a11, a12, a13, a14, a21, a22, a23, a24, a31, ...}
///
void matmult(float* A,float* B,float *C);

void warpShort(short* segw,short* seg1,float* X,int m,int n,int o);

void warpAffine(float* warped,float* input,float* X,int m,int n,int o);

void estimateAffine2(float* X,float* Xprev,float* im1,float* im2,float* costall,float* costall2,int step1,float quant1,int hw1);

#endif // AFFINE_LT1_H_INCLUDED
