#include "dataCostD.h"
#include "transformations.h"
#include <math.h>
#include <algorithm>
#include <iostream>

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

void interp3xyz(float* datai,float* data,float* datax,float* datay,int len1,int len2){
    //x-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            int j2=(j+1)/2;
            if(j%2==1){
                for(int i=0;i<len1;i++){
                    datax[i+j*len1+k*len1*len2]=data[i+j2*len1+k*len1*len1];
                }
            }
            else
            for(int i=0;i<len1;i++){
                datax[i+j*len1+k*len1*len2]=0.5*(data[i+j2*len1+k*len1*len1]+data[i+(j2+1)*len1+k*len1*len1]);
            }
        }
    }
    
    
    //y-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            for(int i=0;i<len2;i++){
                int i2=(i+1)/2;
                if(i%2==1)
                datay[i+j*len2+k*len2*len2]=datax[i2+j*len1+k*len1*len2];
                else
                datay[i+j*len2+k*len2*len2]=0.5*(datax[i2+j*len1+k*len1*len2]+datax[i2+1+j*len1+k*len1*len2]);
            }
        }
    }
    
    //z-interp
    for(int k=0;k<len2;k++){
        int k2=(k+1)/2;
        if(k%2==1){
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=datay[i+j*len2+k2*len2*len2];
                }
            }
            
        }
        else{
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=0.5*(datay[i+j*len2+k2*len2*len2]+datay[i+j*len2+(k2+1)*len2*len2]);
                }
            }
        }
    }
    
}

void interp3xyzB(float* datai,float* data,float* datax,float* datay,int len1,int len2){
    //x-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            int j2=(j+1)/2;
            if(j%2==0){
                for(int i=0;i<len1;i++){
                    datax[i+j*len1+k*len1*len2]=data[i+j2*len1+k*len1*len1];
                }
            }
            else
            for(int i=0;i<len1;i++){
                datax[i+j*len1+k*len1*len2]=0.5*(data[i+j2*len1+k*len1*len1]+data[i+(j2-1)*len1+k*len1*len1]);
            }
        }
    }
    
    
    //y-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            for(int i=0;i<len2;i++){
                int i2=(i+1)/2;
                if(i%2==0)
                datay[i+j*len2+k*len2*len2]=datax[i2+j*len1+k*len1*len2];
                else
                datay[i+j*len2+k*len2*len2]=0.5*(datax[i2+j*len1+k*len1*len2]+datax[i2-1+j*len1+k*len1*len2]);
            }
        }
    }
    
    //z-interp
    for(int k=0;k<len2;k++){
        int k2=(k+1)/2;
        if(k%2==0){
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=datay[i+j*len2+k2*len2*len2];
                }
            }
            
        }
        else{
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=0.5*(datay[i+j*len2+k2*len2*len2]+datay[i+j*len2+(k2-1)*len2*len2]);
                }
            }
        }
    }
    
}


void dataCostCL(unsigned long* data,unsigned long* data2,float* results,int m,int n,int o,int len2,int step1,int hw,float quant,float alpha,int randnum){
    std::cout<<"d"<<std::flush;
    
    int len=hw*2+1;
    len2=pow(hw*2+1,3);
    
    int sz=m*n*o;
    int m1=m/step1; int n1=n/step1; int o1=o/step1;
    int sz1=m1*n1*o1;
    
    //cout<<"len2: "<<len2<<" sz1= "<<sz1<<"\n";
    
    
    
    int quant2=quant;
    
    //const int hw2=hw*quant2; == pad1
    
    int pad1=quant2*hw; int pad2=pad1*2;
    
    int mp=m+pad2; int np=n+pad2; int op=o+pad2;
    int szp=mp*np*op;
    unsigned long* data2p=new unsigned long[szp];
    
    for(int k=0;k<op;k++){
        for(int j=0;j<np;j++){
            for(int i=0;i<mp;i++){
                data2p[i+j*mp+k*mp*np]=data2[std::max(std::min(i-pad1,m-1),0)+std::max(std::min(j-pad1,n-1),0)*m+std::max(std::min(k-pad1,o-1),0)*m*n];
            }
        }
    }
    
    
    int skipz=1; int skipx=1; int skipy=1;
    if(step1>4){
        if(randnum>0){
            skipz=2; skipx=2;
        }
        if(randnum>1){
            skipy=2;
        }
    }
    if(randnum>1&step1>7){
        skipz=3; skipx=3; skipy=3;
    }
    if(step1==4&randnum>1)
    skipz=2;
    
    
    float maxsamp=ceil((float)step1/(float)skipx)*ceil((float)step1/(float)skipz)*ceil((float)step1/(float)skipy);
    //printf("randnum: %d, maxsamp: %d ",randnum,(int)maxsamp);
    
    
    float alphai=(float)step1/(alpha*(float)quant);
    
    float alpha1=0.5*alphai/(float)(maxsamp);
    
    //unsigned long buffer[1000];
    
#pragma omp parallel for
    for(int z=0;z<o1;z++){
        for(int x=0;x<n1;x++){
            for(int y=0;y<m1;y++){
                int z1=z*step1; int x1=x*step1; int y1=y*step1;
                /*for(int k=0;k<step1;k++){
                    for(int j=0;j<step1;j++){
                        for(int i=0;i<step1;i++){
                            buffer[i+j*step1+k*step1*step1]=data[i+y1+(j+x1)*m+(k+z1)*m*n];
                        }
                    }
                }*/
                
                for(int l=0;l<len2;l++){
                    int out1=0;
                    int zs=l/(len*len); int xs=(l-zs*len*len)/len; int ys=l-zs*len*len-xs*len;
                    zs*=quant; xs*=quant; ys*=quant;
                    int x2=xs+x1; int z2=zs+z1; int y2=ys+y1;
                    for(int k=0;k<step1;k+=skipz){
                        for(int j=0;j<step1;j+=skipx){
                            for(int i=0;i<step1;i+=skipy){
                                //unsigned int t=buffer[i+j*STEP+k*STEP*STEP]^buf2p[i+j*mp+k*mp*np];
                                //out1+=(wordbits[t&0xFFFF]+wordbits[t>>16]);
                                unsigned long t1=data[i+y1+(j+x1)*m+(k+z1)*m*n];//buffer[i+j*step1+k*step1*step1];
                                unsigned long t2=data2p[i+j*mp+k*mp*np+(y2+x2*mp+z2*mp*np)];
                                out1+=__builtin_popcountll(t1^t2);
                            }
                        }
                    }
                    results[(y+x*m1+z*m1*n1)*len2+l]=out1*alpha1;
                    
                }
                
            }
        }
    }
    
    
    delete data2p;
    
    return;
    
    
}


void warpImageCL(float* warped, const float* im1, const float* im1b, const float* u1, const float* v1, const float* w1){
    const int m=image_m;
    const int n=image_n;
    const int o=image_o;
    
    float ssd=0;
    float ssd0=0;
    
    interp3(warped,im1,u1,v1,w1,m,n,o,m,n,o,true);
    
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<o;k++){
                ssd+=pow(im1b[i+j*m+k*m*n]-warped[i+j*m+k*m*n],2);
                ssd0+=pow(im1b[i+j*m+k*m*n]-im1[i+j*m+k*m*n],2);
            }
        }
    }
    
    ssd/=m*n*o;
    ssd0/=m*n*o;
    SSD0=ssd0;
    SSD1=ssd;
    
}

void warpAffineS(short* warped,short* input,float* X,float* u1,float* v1,float* w1){
    int m=image_m;
    int n=image_n;
    int o=image_o;
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3]+v1[i+j*m+k*m*n];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7]+u1[i+j*m+k*m*n];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11]+w1[i+j*m+k*m*n];
                int x=round(x1); int y=round(y1);  int z=round(z1);
                
                //if(y>=0&x>=0&z>=0&y<m&x<n&z<o){
                    warped[i+j*m+k*m*n]=input[std::min(std::max(y,0),m-1)+std::min(std::max(x,0),n-1)*m+std::min(std::max(z,0),o-1)*m*n];
                //}
                //else{
                //    warped[i+j*m+k*m*n]=0;
                //}
            }
        }
    }
    
    
}
void warpAffine(float* warped, const float* input, const float * im1b, const float *X, const float *u1, const float *v1, const float *w1){
    // Accepts a 4x4 affine transformation X as input (or at least the first three rows necessary to
    // calculate the transformation in a computer)
    //
    // Accepts two images (of size image_m, image_n, image_o--which are global variables) as inputs:
    // - the image named "input" will be warped using the affine transform X and the displacements u1, v1, w1
    // - the image named "im1b" will be used to measure the sum of squared differences (SSD) before and
    //   after the transformation
    //
    // Accepts three displacement vectors, each having the same dimensions as the input images. These
    // displacement vectors effectively shift each point in the input image elsewhere in the warped
    // image. These shifts, in general (i.e. unless they are constants across all voxels), mean that
    // we're not actually applying an affine transform to the input image. The software elsewhere
    // refers to these displacement vectors as "flow fields" in analogy to fluid mechanics.
    //
    // Accepts an output parameter named "warped," which will store the result of applying the
    // inverse of X to the image named "input."
    //
    // This function does multiple things. It applies an affine transform to the image "input," and
    // it also records the sum of squared errors between the target image im1b and the warped image,
    // as well as between the target image im1b and the original input image. This allows code
    // outside of this function to report the sum of squared differences before and after the
    // transformation, in order to determine how significantly the transformation reduces the error
    // between the two images' intensities.
    //
    // The name of the function is slightly misleading, because it does not warp the image with the
    // affine transform X. Instead, the function applies the inverse of X to the input in order to
    // produce the warped image.
    //
    // Suggested improvements: Split the function into two separate functions (one for transforming
    // the image and one for calculating the errors). Rename the function to reflect what it does.
    //
    
    int m=image_m;
    int n=image_n;
    int o=image_o;
    
    float ssd=0;
    float ssd0=0;
    
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                const float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3]+v1[i+j*m+k*m*n];
                const float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7]+u1[i+j*m+k*m*n];
                const float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11]+w1[i+j*m+k*m*n];
                const int x=floor(x1); const int y=floor(y1);  const int z=floor(z1);
                const float dx=x1-x; const float dy=y1-y; const float dz=z1-z;
                
                
                warped[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*dz*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*dz*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n];
            }
        }
    }
    
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<o;k++){
                ssd+=pow(im1b[i+j*m+k*m*n]-warped[i+j*m+k*m*n],2);
                ssd0+=pow(im1b[i+j*m+k*m*n]-input[i+j*m+k*m*n],2);
            }
        }
    }
    
    ssd/=m*n*o;
    ssd0/=m*n*o;
    SSD0=ssd0;
    SSD1=ssd;
    
    
}

