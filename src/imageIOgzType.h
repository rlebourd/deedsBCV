/*
  Adapted from Mattias P. Heinrich
  Universitaet Luebeck, 2014
 */

#ifndef IMAGE_IOGZ_TYPE_H_INCLUDED
#define IMAGE_IOGZ_TYPE_H_INCLUDED

#include <string>

void writeNifti(std::string filestr,float* vol,char* header,int m,int n,int o,int k);

void writeSegment(std::string filestr,short* vol,char* header,int m,int n,int o);

void gzWriteNifti(std::string filestr,float* data,char* header,int m,int n,int o,int p);

void gzWriteSegment(std::string filestr,short* data,char* header,int m,int n,int o,int p);

template <typename Type>
void readNifti(std::string filestr,Type*& vol,char*& header,int& m,int& n,int& o,int& p){
    
    //convert input std::string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    gzFile file=gzopen(filename,"rb");
    if(!file){
        printf("gzopen of '%s' failed.\n",filename);
        exit(-1);
    }else{
        //cout<<"reading "<<filestr<<"\n";
    }
    
    gzread(file,header,352);
    

    short* dimensions;

    dimensions=reinterpret_cast<short*>(header+40);
    m=(int)dimensions[1];
    n=(int)dimensions[2];
    o=(int)dimensions[3];
    p=(int)dimensions[4];
    
    int sz=m*n*o*p;
    
    vol=new Type[sz];

    //printf("Read image with dimensions %dx%dx%dx%d\n",m,n,o,p);
    short* datatype; //read datatype (and bitpix)
    datatype=reinterpret_cast<short*>(header+70);
    short bitpix=datatype[1];
    
    long charbytes=bitpix/8;
    //
    char* filecharptr=new char[m*n*o*p*charbytes];
    //raw empty datapointers of all 'types'
    unsigned char* ucharptr; float* floatptr; short* shortptr;
    unsigned short* ushortptr; int* intptr; double* doubleptr;
    //read input data values depending on datatype and convert to float
    //binary are read character by character, reinterpret_cast converts them
    //copy them into float* array afterwards
    switch(datatype[0]){
        case 2:
            gzread(file,filecharptr,m*n*o*p*sizeof(unsigned char));
            ucharptr=reinterpret_cast<unsigned char*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=ucharptr[i];
            }
            break;
        case 4:
            gzread(file,filecharptr,m*n*o*p*sizeof(short));
            shortptr=reinterpret_cast<short*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=shortptr[i];
            }
            break;
        case 8:
            gzread(file,filecharptr,m*n*o*p*sizeof(int));
            intptr=reinterpret_cast<int*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=intptr[i];
            }
            break;
        case 16:
            gzread(file,filecharptr,m*n*o*p*sizeof(float));
            floatptr=reinterpret_cast<float*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=floatptr[i];
            }
            break;
        case 64:
            gzread(file,filecharptr,m*n*o*p*sizeof(double));
            doubleptr=reinterpret_cast<double*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=doubleptr[i];
            }
            break;
        case 512:
            gzread(file,filecharptr,m*n*o*p*sizeof(unsigned short));
            ushortptr=reinterpret_cast<unsigned short*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=ushortptr[i];
            }
            break;
        default:
            printf("Datatype %d not supported. Exiting.\n",datatype[0]);
            exit(1);
    }
    delete filecharptr;
    
    //gzread(file,vol,sz*sizeof(Type));
    //vol=reinterpret_cast<Type*>(vol);
    
    gzclose(file);
    
}

template <typename TypeFO>
void writeOutput(TypeFO* data,std::string filestr,int length){
    //opens file for binary-output
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    ofstream ofs1(filename,ios::out|ios::binary);
    ofs1.write(reinterpret_cast<char*>(data),length*sizeof(TypeFO));
    ofs1.close();
}

template <typename TypeF>
vector<TypeF> readFile(std::string filestr){
    //convert input std::string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    //opens file for binary-input
    ifstream file(filename,ios::binary|ios::ate);
    vector<TypeF> invalues;
    if(file.is_open()){
        int length=file.tellg()/(sizeof(TypeF));
        file.seekg(0,file.beg); //go to beginning
        char* charptr=new char[length*sizeof(TypeF)];
        file.read((char*)charptr,length*sizeof(TypeF));
        TypeF* floatptr=reinterpret_cast<TypeF*>(charptr);
        for(int i=0;i<length;i++){
            invalues.push_back(floatptr[i]);
        }
        delete charptr;
    }
    else{
        cout<<"readFile error. Did not find "<<filestr<<"\n";
        exit(-1);
    }
    file.close();
    return invalues;
}

#endif // IMAGE_IOGZ_TYPE_H_INCLUDED
