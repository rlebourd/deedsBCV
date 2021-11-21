#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include <string.h>
#include <sstream>
#include <pthread.h>
#include <thread>
#include "zlib.h"
#include <sys/stat.h>

using namespace std;
//compile with openMP g++ linearBCV.cpp -O3 -std=c++11 -mavx2 -msse4.2 -pthread -fopenmp -lz -o linear11

//some global variables
int RAND_SAMPLES; //will all be set later (if needed)
int image_m; int image_n; int image_o; int image_d=1;
float SSD0=0.0; float SSD1=0.0; float SSD2=0.0; float distfx_global; float beta=1;
//float SIGMA=8.0;
int qc=1;

//struct for multi-threading of mind-calculation
struct mind_data{
	float* im1;
    float* d1;
    uint64_t* mindq;
    int qs;
    int ind_d1;
};

float quantile(float* values,int length,float quant1){
    float* values2=new float[length];
    for(int i=0;i<length;i++){
        values2[i]=values[i];
    }
    int quantind=length*min(max(quant1,0.01f),0.99f);
    //printf("quantind: %d/%d\n",quantind,length);
    nth_element(values2,values2+quantind,values2+length);
    float med1=values2[quantind];
    // delete values2;
    return med1;
}
bool RIGID=false;

#include "imageIOgzType.h"
#include "transformations.h"
#include "QRsolve.h"
#include "affineLTS1.h"
#include "MINDSSCbox.h"
#include "dataCostD.h"
#include "parseArguments.h"

static void printUsage(){
    cout<<"=============================================================\n";
    cout<<"Usage (required input arguments):\n";
    cout<<"./linearBCV -F fixed.nii.gz -M moving.nii.gz -O output\n";
    cout<<"optional parameters:\n";
    cout<<" -R <find rigid instead of affine transform> (default 0)\n";
    cout<<" -l <number of levels> (default 4)\n";
    cout<<" -G <grid spacing for each level> (default 7x6x5x4)\n";
    cout<<" -L <maximum search radius - each level> (default 5x4x3x2)\n";
    cout<<" -Q <quantisation of search step size> (default 4x3x2x1)\n";
    cout<<"=============================================================\n";
}

static bool shouldPrintUsageBasedOnArgs(int argc, char * const argv[]){
    return (argc<4||argv[1][1]=='h');
}

static bool fileHasCorrectExtension(std::string filename){
    return filename.size() >= 2 && filename.substr(filename.size() - 2) == "gz";
}

static std::string getLastPathComponent(std::string filepath){
    std::string forewardOrBackwardSlash{"/\\"};
    const auto iter = filepath.find_last_of(forewardOrBackwardSlash);
    if (iter != filepath.npos){
        return filepath.substr(iter+1);
    } else {
        return filepath;
    }
}

int main (int argc, char * const argv[]) {
    //PARSE INPUT ARGUMENTS
    
    // 1. validate inputs
    if(shouldPrintUsageBasedOnArgs(argc, argv)){
        float Xinv[16]; float Ident[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
        float Xprev[16]={1,2,3,4,0,1,10,100,1,2,1,6,7,10,0,1};
        qrsolve(Xinv,Xprev,Ident,4,4);
        matmult(Xinv, Xprev, Ident);
        printUsage();
        return 1;
    }
    
    // 2. parse program's inputs
    parameters args{
        //defaults
        .grid_spacing{7,6,5,4},
        .search_radius{5,4,3,2},
        .quantisation{4,3,2,1},
        .levels = 4,
    };
    parseCommandLine(args, argc, argv);
    
    if(fileHasCorrectExtension(args.fixed_file) != true){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    if(fileHasCorrectExtension(args.moving_file) != true){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    
    // 3. print out that we're beginning registration
    cout<<"Starting linear registration of "<<getLastPathComponent(args.fixed_file)<<" and "<<getLastPathComponent(args.moving_file)<<"\n";
    cout<<"=============================================================\n";

    // 4. various, apparently unrelated parameters are set
    if(args.rigid){
        cout<<"Searching only for rigid (6 parameter) transform\n";
        RIGID=true;
    }

	qc=2;
    
	float alpha=1;
    //READ IMAGES and INITIALISE ARRAYS
        
    RAND_SAMPLES=1; //fixed/efficient random sampling strategy
    
    // 5. read in the images from the specified files and validate that their sizes are compatible
    float* movingImageData; float* fixedImageData;
    
    int M,N,O,P; //image dimensions
    
    //==ALWAYS ALLOCATE MEMORY FOR HEADER ===/
    char* header=new char[352];
    
    readNifti(args.fixed_file,fixedImageData,header,M,N,O,P);
    image_m=M; image_n=N; image_o=O;

    readNifti(args.moving_file,movingImageData,header,M,N,O,P);
    
    if (M!=image_m|N!=image_n|O!=image_o){
        cout<<"Inconsistent image sizes (must have same dimensions) " << M << " " << N << " " << " " << O << " and " << image_m << " " << image_n << " " << image_o << endl;
        return -1;
    }
    cout<<"Consistent image sizes (must have same dimensions) " << M << " " << N << " " << " " << O << endl;
    
    // 6. scale the image's intensity values by 1024 (for some reason this is needed for CTs)
    int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
    
    //assume we are working with CT scans (add 1024 HU)
    float thresholdF=-1024; float thresholdM=-1024;
    
    for (int i=0;i<sz;i++){
        fixedImageData[i]-=thresholdF;
        movingImageData[i]-=thresholdM;
    }
    
    // 7.
    float *movingImageDataWarped=new float[m*n*o];
    float *fixedImageDataWarped=new float[m*n*o];

    vector<int> mind_step;
    for (int i=0;i<args.quantisation.size();i++){
        mind_step.push_back(floor(0.5f*(float)args.quantisation[i]+1.0f));
        // for quantisations {4, 3, 2, 1}, the mind_steps become {3, 2, 2, 1}
    }
    float* X=new float[4*4];
    
    float Xprev[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    for (int i=0;i<16;i++){
        X[i]=Xprev[i];
    }
    
    uint64_t *movingImageMindDescData = new uint64_t[m*n*o];
    uint64_t *fixedImageMindDescData = new uint64_t[m*n*o];
    uint64_t *warpedImageMindDescData = new uint64_t[m*n*o];

    //==========================================================================================
    //==========================================================================================
    const float Ident[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    float Xinv[16];
    for (int level=0; level<args.levels; level++){
        // initializes the parameters for the current iteration
        const float quantisation = args.quantisation[level];
        const int gridSpacing = args.grid_spacing[level];
        const int searchRadius = args.search_radius[level];

        // finds the matrix inverse of Xprev
        qrsolve(Xinv,Xprev,Ident,4,4);
        
        //
        warpAffine(fixedImageDataWarped,fixedImageData,Xinv,m,n,o); // so now fixedImageDataWarped should be close to movingImageData
        warpAffine(movingImageDataWarped,movingImageData,Xprev,m,n,o); // so now movingImageDataWarped should be close to fixedImageData

        // update descriptors
        const float prevMindStep = mind_step[max(level-1, 0)];
        const float currMindStep = mind_step[level];
        
        if(level == 0 || prevMindStep != currMindStep){
            // the descriptors are always calculated on the first iteration
            // and then re-calculated whenever the "mind step" changes
            
            descriptor(movingImageMindDescData,movingImageData,m,n,o,currMindStep); // so now movingImageMindDescData contains descriptor data for movingImageData
            descriptor(fixedImageMindDescData,fixedImageData,m,n,o,currMindStep); // so now fixedImageMindDescData contains descriptor data for fixedImageData
        }

        // forward and reverse registration
        const int len3 = pow(searchRadius*2+1, 3); // (side length of each patch)^3 = volume of patch
        const int m1 = m/gridSpacing; const int n1 = n/gridSpacing; const int o1 = o/gridSpacing; const int sz1 = m1*n1*o1;

        //FULL-REGISTRATION FORWARDS
        float* costall = new float[sz1*len3];
        descriptor(warpedImageMindDescData,movingImageDataWarped,m,n,o,currMindStep); // so now warpedImageMindDescData should contain descriptors that resemble those of the fixedImageMindDescData
        dataCostCL((unsigned long*)fixedImageMindDescData,(unsigned long*)warpedImageMindDescData,costall,m,n,o,len3,gridSpacing,searchRadius,quantisation,alpha,RAND_SAMPLES);

        //FULL-REGISTRATION BACKWARDS
        float* costall2 = new float[sz1*len3];
        descriptor(warpedImageMindDescData,fixedImageDataWarped,m,n,o,currMindStep); // so now warpedImageMindDescData should contain descriptors that resemble those of the movingImageMindDescData
        dataCostCL((unsigned long*)movingImageMindDescData,(unsigned long*)warpedImageMindDescData,costall2,m,n,o,len3,gridSpacing,searchRadius,quantisation,alpha,RAND_SAMPLES);
        
        estimateAffine2(X,Xprev,fixedImageData,movingImageData,costall,costall2,gridSpacing,quantisation,searchRadius);
        
        delete[] costall;
        delete[] costall2;

        // update the affine transformation matrix
        for (int i=0; i<16; i++){
            Xprev[i] = X[i];
        }
	}
    delete[] movingImageMindDescData;
    delete[] fixedImageMindDescData;
    
    //==========================================================================================
    //==========================================================================================
    
    string outputfile;
    outputfile.append(args.output_stem);
    outputfile.append("_matrix.txt");

    ofstream matfile;
    matfile.open(outputfile);
    
    for (int i=0;i<4;i++){
        matfile<<X[i]<<"  "<<X[i+4]<<"  "<<X[i+8]<<"  "<<X[i+12]<<"\n";
        
        printf("%+4.3f | %+4.3f | %+4.3f | %+4.3f \n",X[i],X[i+4],X[i+8],X[i+12]);
    }
    matfile.close();
    
    // if SEGMENTATION of moving image is provided APPLY SAME TRANSFORM
    if (args.segment){
        short* seg2;
        readNifti(args.moving_seg_file,seg2,header,M,N,O,P);
        
        float* zero=new float[sz];
        fill(zero,zero+sz,0.0f);

        short* segw=new short[sz];
        fill(segw,segw+sz,0);
        
        warpAffineS(segw,seg2,X,zero,zero,zero);
        
        string outputseg;
        outputseg.append(args.output_stem);
        outputseg.append("_deformed_seg.nii.gz");

        gzWriteSegment(outputseg,segw,header,m,n,o,1);
    }
  
    return 0;
}
