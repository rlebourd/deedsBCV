#ifndef QR_SOLVE_H_INCLUDED
#define QR_SOLVE_H_INCLUDED

float norm(float* vector,int len);

float dotprod(float* vector,float* vector2,int len);

void qrsolve(float* X,float* A,float* b,int len,int len2);

void jacobiSVD3(float* A,float* U,float* V);

void findRigid(float* RT,float* pts1t,float* pts2t,int len);

void affineRobust(float* RT,float* pts1,float* pts2,int len);

#endif // QR_SOLVE_H_INCLUDED
