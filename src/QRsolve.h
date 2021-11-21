#ifndef QR_SOLVE_H_INCLUDED
#define QR_SOLVE_H_INCLUDED

float norm(float* vector,int len);

float dotprod(const float* vector,const float* vector2,int len);

/// @brief Finds the solution to a system of linear equations
///
/// It's not abundantly clear whether the len and len2 arguments actually matter.
/// The program passes them in both as the number 4, and the algorithm as
/// written appears to have been written assuming a 4x4 matrix, so I think these
/// parameters are essentially a false promise in the interface.
///
/// @param [out] X is the matrix inverse of the 4x4 matrix A
/// @param [in] A is the 4x4 matrix to invert
/// @param [in] b is the identity matrix in the examples
///
void qrsolve(float* X,float* A,const float* b,int len,int len2);

void jacobiSVD3(float* A,float* U,float* V);

void findRigid(float* RT,float* pts1t,float* pts2t,int len);

void affineRobust(float* RT,float* pts1,float* pts2,int len);

#endif // QR_SOLVE_H_INCLUDED
