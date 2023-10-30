#include <fstream>
#include <stdio.h>
#include <stdexcept>
#include <cstdlib>
#include <math.h>
#include <algorithm>

#include "float2DReg.h"
#include "doubleHyperExt.h"
#include "axis.h"
#include "hypercube.h"
#include "seplib.h"


using namespace SEP;

void VECDEXT::sepWrite(const std::shared_ptr<doubleHyper> vec, std::string output){
    unsigned long long n123 = vec->getHyper()->getN123();
    std::vector<axis> axes = vec->getHyper()->getAxes();
    std::string ni, oi, di, labeli;
    const char* out=output.c_str();
    for (int i=0; i<axes.size(); i++){
        ni = "n" + std::to_string(i+1);
        oi = "o" + std::to_string(i+1);
        di = "d" + std::to_string(i+1);
        labeli = "label" + std::to_string(i+1);

        if(0 != auxputch(ni.c_str(),"d",&axes[i].n,out)){
            ni = "Error: Cannot write n" + std::to_string(i+1)+"\n";
    	    seperr(ni.c_str());
        }

        if(0 != auxputch(oi.c_str(),"f",&axes[i].o,out)){
            oi = "Error: Cannot write o" + std::to_string(i+1)+"\n";
            seperr(oi.c_str());
        }
    	    

        if(0 != auxputch(di.c_str(),"f",&axes[i].d,out)){
            di = "Error: Cannot write d" + std::to_string(i+1)+"\n";
            seperr(di.c_str());
        }   

        if(0 != auxputch(labeli.c_str(),"s",axes[i].label.c_str(),out)){
            labeli = "Error: Cannot write label" + std::to_string(i+1)+"\n";
    	    seperr(labeli.c_str());
        }
    }
	
	float* val = new float[n123];
	for (int i=0; i<n123; i++){
		val[i] = vec->getVals()[i];
	}

	if(4*n123 != srite(out,val,4*n123))
		seperr("Error: Cannot write data\n");
	
	delete val;
}


std::shared_ptr<double1DReg> VECDEXT::sepRead1D(){
    int n1;
	float o1, d1;
	char label1[20];
	if(0 == hetch((char*)"n1",(char*)"d",&n1))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == hetch((char*)"o1",(char*)"f",&o1))
    	seperr("Error: Cannot read o1\n");

	if(0 == hetch((char*)"d1",(char*)"f",&d1))
    	seperr("Error: Cannot read d1\n");

	if(0 == hetch((char*)"label1",(char*)"s",&label1))
    	seperr("Error: Cannot read label1\n");
	
	std::string lab1(label1);

	axis X(n1, o1, d1, lab1);
	unsigned long long n123 = X.n;
	std::shared_ptr<double1DReg> vec (new double1DReg(X));

	float* val = new float[n123];

	if(4*n123 != sreed("in",val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

std::shared_ptr<double1DReg> VECDEXT::sepRead1D(std::string input){
    int n1;
	float o1, d1;
	char label1[20];
    const char* data = input.c_str();
	if(0 == auxpar((char*)"n1",(char*)"d",&n1,data))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == auxpar((char*)"o1",(char*)"f",&o1,data))
    	seperr("Error: Cannot read o1\n");

	if(0 == auxpar((char*)"d1",(char*)"f",&d1,data))
    	seperr("Error: Cannot read d1\n");

	if(0 == auxpar((char*)"label1",(char*)"s",&label1,data))
    	seperr("Error: Cannot read label1\n");
	
	std::string lab1(label1);

	axis Z(n1, o1, d1, lab1);
	unsigned long long n123 = Z.n;
	std::shared_ptr<double1DReg> vec (new double1DReg(Z));

	float* val = new float[n123];

	if(4*n123 != sreed(data,val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

std::shared_ptr<double2DReg> VECDEXT::sepRead(){
    int n1, n2;
	float o1, o2, d1, d2;
	char label1[20]; char label2[20];
	if(0 == hetch((char*)"n1",(char*)"d",&n1))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == hetch((char*)"n2",(char*)"d",&n2))
    	seperr("Error: Cannot read n2\n");
	
	if(0 == hetch((char*)"o1",(char*)"f",&o1))
    	seperr("Error: Cannot read o1\n");

	if(0 == hetch((char*)"o2",(char*)"f",&o2))
    	seperr("Error: Cannot read o2\n");

	if(0 == hetch((char*)"d1",(char*)"f",&d1))
    	seperr("Error: Cannot read d1\n");

	if(0 == hetch((char*)"d2",(char*)"f",&d2))
    	seperr("Error: Cannot read d2\n");

	if(0 == hetch((char*)"label1",(char*)"s",&label1))
    	seperr("Error: Cannot read label1\n");

	if(0 == hetch((char*)"label2",(char*)"s",&label2))
    	seperr("Error: Cannot read label2\n");
	
	std::string lab1(label1);
	std::string lab2(label2);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
	unsigned long long n123 = X.n*Z.n;
	std::shared_ptr<double2DReg> vec (new double2DReg(Z, X));

	float* val = new float[n123];

	if(4*n123 != sreed("in",val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

std::shared_ptr<double2DReg> VECDEXT::sepRead(std::string input){
    int n1, n2;
	float o1, o2, d1, d2;
	char label1[20]; char label2[20];
    const char* data = input.c_str();
	if(0 == auxpar((char*)"n1",(char*)"d",&n1,data))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == auxpar((char*)"n2",(char*)"d",&n2,data))
    	seperr("Error: Cannot read n2\n");
	
	if(0 == auxpar((char*)"o1",(char*)"f",&o1,data))
    	seperr("Error: Cannot read o1\n");

	if(0 == auxpar((char*)"o2",(char*)"f",&o2,data))
    	seperr("Error: Cannot read o2\n");

	if(0 == auxpar((char*)"d1",(char*)"f",&d1,data))
    	seperr("Error: Cannot read d1\n");

	if(0 == auxpar((char*)"d2",(char*)"f",&d2,data))
    	seperr("Error: Cannot read d2\n");

	if(0 == auxpar((char*)"label1",(char*)"s",&label1,data))
    	seperr("Error: Cannot read label1\n");

	if(0 == auxpar((char*)"label2",(char*)"s",&label2,data))
    	seperr("Error: Cannot read label2\n");
	
	std::string lab1(label1);
	std::string lab2(label2);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
	unsigned long long n123 = X.n*Z.n;

	float * val = new float[n123];
	std::shared_ptr<double2DReg> vec (new double2DReg(Z, X));

	if(4*n123 != sreed(data,val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

std::shared_ptr<double3DReg> VECDEXT::sepRead3D(){
    int n1, n2, n3;
	float o1, o2, o3, d1, d2, d3;
	char label1[20]; char label2[20]; char label3[20];
	if(0 == hetch((char*)"n1",(char*)"d",&n1))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == hetch((char*)"n2",(char*)"d",&n2))
    	seperr("Error: Cannot read n2\n");

    if(0 == hetch((char*)"n3",(char*)"d",&n3))
    	seperr("Error: Cannot read n3\n");
	
	if(0 == hetch((char*)"o1",(char*)"f",&o1))
    	seperr("Error: Cannot read o1\n");

	if(0 == hetch((char*)"o2",(char*)"f",&o2))
    	seperr("Error: Cannot read o2\n");

    if(0 == hetch((char*)"o3",(char*)"f",&o3))
    	seperr("Error: Cannot read o3\n");

	if(0 == hetch((char*)"d1",(char*)"f",&d1))
    	seperr("Error: Cannot read d1\n");

	if(0 == hetch((char*)"d2",(char*)"f",&d2))
    	seperr("Error: Cannot read d2\n");

    if(0 == hetch((char*)"d3",(char*)"f",&d3))
    	seperr("Error: Cannot read d3\n");

	if(0 == hetch((char*)"label1",(char*)"s",&label1))
    	seperr("Error: Cannot read label1\n");

	if(0 == hetch((char*)"label2",(char*)"s",&label2))
    	seperr("Error: Cannot read label2\n");

    if(0 == hetch((char*)"label3",(char*)"s",&label3))
    	seperr("Error: Cannot read label3\n");
	
	std::string lab1(label1);
	std::string lab2(label2);
    std::string lab3(label3);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
    axis Y(n3, o3, d3, lab3);
	unsigned long long n123 = Y.n*X.n*Z.n;
	std::shared_ptr<double3DReg> vec (new double3DReg(Z, X, Y));

	float* val = new float[n123];

	if(4*n123 != sreed("in",val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

std::shared_ptr<double3DReg> VECDEXT::sepRead3D(std::string input){
    int n1, n2, n3;
	float o1, o2, o3, d1, d2, d3;
	char label1[20]; char label2[20]; char label3[20];
    const char* data = input.c_str();
	if(0 == auxpar((char*)"n1",(char*)"d",&n1,data))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == auxpar((char*)"n2",(char*)"d",&n2,data))
    	seperr("Error: Cannot read n2\n");

    if(0 == auxpar((char*)"n3",(char*)"d",&n3,data))
    	seperr("Error: Cannot read n3\n");
	
	if(0 == auxpar((char*)"o1",(char*)"f",&o1,data))
    	seperr("Error: Cannot read o1\n");

	if(0 == auxpar((char*)"o2",(char*)"f",&o2,data))
    	seperr("Error: Cannot read o2\n");

    if(0 == auxpar((char*)"o3",(char*)"f",&o3,data))
    	seperr("Error: Cannot read o3\n");

	if(0 == auxpar((char*)"d1",(char*)"f",&d1,data))
    	seperr("Error: Cannot read d1\n");

	if(0 == auxpar((char*)"d2",(char*)"f",&d2,data))
    	seperr("Error: Cannot read d2\n");

    if(0 == auxpar((char*)"d3",(char*)"f",&d3,data))
    	seperr("Error: Cannot read d3\n");

	if(0 == auxpar((char*)"label1",(char*)"s",&label1,data))
    	seperr("Error: Cannot read label1\n");

	if(0 == auxpar((char*)"label2",(char*)"s",&label2,data))
    	seperr("Error: Cannot read label2\n");

    if(0 == auxpar((char*)"label3",(char*)"s",&label3,data))
    	seperr("Error: Cannot read label3\n");
	
	std::string lab1(label1);
	std::string lab2(label2);
    std::string lab3(label3);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
    axis Y(n3, o3, d3, lab3);
	unsigned long long n123 = Y.n*X.n*Z.n;
	std::shared_ptr<double3DReg> vec (new double3DReg(Z, X, Y));

	float* val = new float[n123];

	if(4*n123 != sreed(data,val,4*n123))
    	seperr("Error: Cannot read data\n");

	for (unsigned long i=0; i<n123; i++){
		vec->getVals()[i] = val[i];
	}

	delete val;

    return vec;
}

void VECDEXT::random(std::shared_ptr<doubleHyper> vec, float min, float max){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();
    for (unsigned long long i=0; i<n123; i++){
            vals[i] += min + (max-min)*(rand() % 10000)/10000.0;
    }    
}

double VECDEXT::norm2(const std::shared_ptr<doubleHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i]*vals[i];
    } 
    return sum;   
}

double VECDEXT::sum(const std::shared_ptr<doubleHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i];
    } 
    return sum;   
}

double VECDEXT::rms(const std::shared_ptr<doubleHyper> vec){
    return sqrt(norm2(vec)/(vec->getHyper()->getN123()));
}