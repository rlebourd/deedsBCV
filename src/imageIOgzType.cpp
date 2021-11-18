#include "imageIOgzType.h"
#include <string>

void writeNifti(string filestr,float* vol,char* header,int m,int n,int o,int k){
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(k>1){
        dim=4;
    }
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n; dimensions[3]=o;
    dimensions[4]=k;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+8,header2+40);

    short* datatype=new short[2];
    datatype[0]=16; datatype[1]=32; //float datatype
    char* datachar=reinterpret_cast<char*>(datatype);
    copy(datachar,datachar+4,header2+70);
    
    //printf("Writing image with dimensions %dx%dx%d\n",m,n,o);

    ofstream file(filename,ios::out|ios::binary);
    //opens file for binary-output
	if(file.is_open()){
		file.write(header2,352);
		file.write(reinterpret_cast<char*>(vol),m*n*o*4);
		file.close();
		cout<<"File "<<filename<<" written.\n";
	}
    else{
        printf("File error. Could not write file.\n");
    }
}

void writeSegment(string filestr,short* vol,char* header,int m,int n,int o){
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n; dimensions[3]=o;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+8,header2+40);
    
    short* datatype=new short[2];
    datatype[0]=4; datatype[1]=512; //short datatype
    char* datachar=reinterpret_cast<char*>(datatype);
    copy(datachar,datachar+4,header2+70);
    
    //printf("Writing segmentation with dimensions %dx%dx%d\n",m,n,o);
    
    ofstream file(filename,ios::out|ios::binary);
    //opens file for binary-output
	if(file.is_open()){
		file.write(header2,352);
		file.write(reinterpret_cast<char*>(vol),m*n*o*2);
		file.close();
		cout<<"File "<<filename<<" written.\n";
	}
    else{
        printf("File error. Could not write file.\n");
    }
}


void gzWriteNifti(string filestr,float* data,char* header,int m,int n,int o,int p){
    
    int sz=m*n*o*p;
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    //do not change original header, define (new )dimensions and datatype
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(p>1){ dim=4;}
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n;
    dimensions[3]=o; dimensions[4]=p;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+10,header2+40);
    short* datatype=new short[2];
    datatype[0]=16; datatype[1]=32; //float datatype
    char* typechar=reinterpret_cast<char*>(datatype);
    copy(typechar,typechar+4,header2+70);
    
    //copy header and data into out
    int size=352+sz*sizeof(float);
    char* out2=new char[size];
    copy(header2,header2+352,out2);
    char* datachar=reinterpret_cast<char*>(data);
    copy(datachar,datachar+sz*sizeof(float),out2+352);
    
    //write gzipped file
    gzFile file2=gzopen(filename,"wb");
    if(!file2){
        printf("File error. Could not write file.\n");
        exit(1);
    }
    int err=gzwrite(file2,(unsigned char*)out2,size);
    gzclose(file2);
    
    delete out2; delete header2;
}

void gzWriteSegment(string filestr,short* data,char* header,int m,int n,int o,int p){
    int sz=m*n*o*p;
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    //do not change original header, define (new )dimensions and datatype
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(p>1){ dim=4;}
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n;
    dimensions[3]=o; dimensions[4]=p;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+10,header2+40);
    short* datatype=new short[2];
    datatype[0]=4; datatype[1]=16;// edit 08-16 - 512; //short datatype
    char* typechar=reinterpret_cast<char*>(datatype);
    copy(typechar,typechar+4,header2+70);
    
    //copy header and data into out
    int size=352+sz*sizeof(short);
    char* out2=new char[size];
    copy(header2,header2+352,out2);
    char* datachar=reinterpret_cast<char*>(data);
    copy(datachar,datachar+sz*sizeof(short),out2+352);
    
    //write gzipped file
    gzFile file2=gzopen(filename,"wb");
    if(!file2){
        printf("File error. Could not write file.\n");
        exit(1);
    }
    int err=gzwrite(file2,(unsigned char*)out2,size);
    gzclose(file2);
    
    delete out2; delete header2;
}
