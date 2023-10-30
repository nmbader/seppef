#include <fstream>
#include <stdio.h>
#include "floatHyperExt.h"
#include "complexHyperExt.h"
#include "complex2DReg.h"
#include "axis.h"
#include <stdexcept>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include "functions.h"
#include "hypercube.h"
#include "axis.h"
#include "fkTransform.h"
#include "fxTransform.h"
#include "seplib.h"


using namespace SEP;

float VECEXT::getox(const std::shared_ptr<float2DReg> vec) {
    return vec->getHyper()->getAxis(2).o;}

float VECEXT::getoz(const std::shared_ptr<float2DReg> vec) {
return vec->getHyper()->getAxis(1).o;}

float VECEXT::getdx(const std::shared_ptr<float2DReg> vec) {
return vec->getHyper()->getAxis(2).d;}

float VECEXT::getdz(const std::shared_ptr<float2DReg> vec) {
return vec->getHyper()->getAxis(1).d;}

unsigned int VECEXT::getnx(const std::shared_ptr<floatHyper> vec) {
return vec->getHyper()->getAxis(2).n;}

unsigned int VECEXT::getnz(const std::shared_ptr<floatHyper> vec) {
return vec->getHyper()->getAxis(1).n;}

void VECEXT::setox(std::shared_ptr<floatHyper> vec, float ox) {
    std::vector<axis> axes = vec->getHyper()->getAxes();
    axes[1].o = ox;
    std::shared_ptr<hypercube> h (new hypercube(axes));
    vec->setHyper(h);
}

void VECEXT::setoz(std::shared_ptr<floatHyper> vec, float oz) {
    std::vector<axis> axes = vec->getHyper()->getAxes();
    axes[0].o = oz;
    std::shared_ptr<hypercube> h (new hypercube(axes));
    vec->setHyper(h);
}

void VECEXT::write(const std::shared_ptr<floatHyper> vec, std::string name) {
    unsigned long long n123 = vec->getHyper()->getN123();
    std::ofstream output;
    output.open(name+".bin", std::ios::out | std::ios::binary);
    if (output.is_open())
	    output.write((char*)vec->getVals(), n123*sizeof(float));
	output.close();

    std::vector<axis> axes = vec->getHyper()->getAxes();
    output.open(name+".header");
    if (output.is_open()){
        output << "type=float"<< axes.size() << "DReg\n";
        for (int i=0; i<axes.size(); i++){
            output << "n"<< i+1 <<"=" << axes[i].n << "\n";
            output << "o"<< i+1 <<"=" << axes[i].o << "\n";
            output << "d"<< i+1 <<"=" << axes[i].d << "\n";
        }
    }
    output.close();
}

std::shared_ptr<hypercube> VECEXT::readHyper(std::string name){
    std::vector<axis> axes;
    int ndim=0;
    std::string line;
    std::ifstream input;
	input.open(name+".header");
    if (input.is_open()){
        getline(input,line);
		line.erase(0,10);
		line.erase(1,4);
        ndim = std::stoi(line.c_str());
    }

    for (int i=0; i<ndim; i++){
        axis ax;
        getline(input,line);
        line.erase(0,3);
        ax.n = std::stoi(line.c_str());
        getline(input,line);
        line.erase(0,3);
        ax.o = std::stof(line.c_str());
        getline(input,line);
        line.erase(0,3);
        ax.d = std::stof(line.c_str());
        axes.push_back(ax);
    }
    input.close();
    std::shared_ptr<hypercube> hyper( new hypercube(axes));
    return hyper;
}

void VECEXT::read(const std::shared_ptr<floatHyper> vec, std::string name){
    unsigned long long n123 = vec->getHyper()->getN123();
    std::ifstream input;
	input.open(name+".bin", std::ios::in | std::ios::binary);
    if (input.is_open())
	    input.read((char*) vec->getVals(), n123*sizeof(float));
	input.close();
}

void VECEXT::sepWrite(const std::shared_ptr<floatHyper> vec, std::string output){
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
    
        if(4*n123 != srite(out,vec->getVals(),4*n123))
    	    seperr("Error: Cannot write data\n");
}

void VECEXT::sepWrite(const std::shared_ptr<doubleHyper> vec, std::string output){
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

std::shared_ptr<float1DReg> VECEXT::sepRead1D(){
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
	std::shared_ptr<float1DReg> vec (new float1DReg(X));

	if(4*n123 != sreed("in",vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

std::shared_ptr<float1DReg> VECEXT::sepRead1D(std::string input){
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
	std::shared_ptr<float1DReg> vec (new float1DReg(Z));

	if(4*n123 != sreed(data,vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

std::shared_ptr<float2DReg> VECEXT::sepRead(){
    int n1, n2;
	float o1, o2, d1, d2;
	char label1[20]; char label2[20];
	if(0 == hetch((char*)"n1",(char*)"d",&n1))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == hetch((char*)"n2",(char*)"d",&n2))
    	//seperr("Error: Cannot read n2\n");
        n2=1;
	
	if(0 == hetch((char*)"o1",(char*)"f",&o1))
    	seperr("Error: Cannot read o1\n");

	if(0 == hetch((char*)"o2",(char*)"f",&o2))
    	//seperr("Error: Cannot read o2\n");
        o2=0;

	if(0 == hetch((char*)"d1",(char*)"f",&d1))
    	seperr("Error: Cannot read d1\n");

	if(0 == hetch((char*)"d2",(char*)"f",&d2))
    	//seperr("Error: Cannot read d2\n");
        d2=1;

	if(0 == hetch((char*)"label1",(char*)"s",&label1))
    	//seperr("Error: Cannot read label1\n");
        {};

	if(0 == hetch((char*)"label2",(char*)"s",&label2))
    	//seperr("Error: Cannot read label2\n");
        {};
	
	std::string lab1(label1);
	std::string lab2(label2);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
	unsigned long long n123 = X.n*Z.n;
	std::shared_ptr<float2DReg> vec (new float2DReg(Z, X));

	if(4*n123 != sreed("in",vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

std::shared_ptr<float2DReg> VECEXT::sepRead(std::string input){
    int n1, n2;
	float o1, o2, d1, d2;
	char label1[20]; char label2[20];
    const char* data = input.c_str();
	if(0 == auxpar((char*)"n1",(char*)"d",&n1,data))
    	seperr("Error: Cannot read n1\n");
	
	if(0 == auxpar((char*)"n2",(char*)"d",&n2,data))
    	//seperr("Error: Cannot read n2\n");
        n2=1;
	
	if(0 == auxpar((char*)"o1",(char*)"f",&o1,data))
    	seperr("Error: Cannot read o1\n");

	if(0 == auxpar((char*)"o2",(char*)"f",&o2,data))
    	//seperr("Error: Cannot read o2\n");
        o2=0;

	if(0 == auxpar((char*)"d1",(char*)"f",&d1,data))
    	seperr("Error: Cannot read d1\n");

	if(0 == auxpar((char*)"d2",(char*)"f",&d2,data))
    	//seperr("Error: Cannot read d2\n");
        d2=1;

	if(0 == auxpar((char*)"label1",(char*)"s",&label1,data))
    	//seperr("Error: Cannot read label1\n");
        {};

	if(0 == auxpar((char*)"label2",(char*)"s",&label2,data))
    	//seperr("Error: Cannot read label2\n");
        {};
	
	std::string lab1(label1);
	std::string lab2(label2);

	axis Z(n1, o1, d1, lab1);
	axis X(n2, o2, d2, lab2);
	unsigned long long n123 = X.n*Z.n;
	std::shared_ptr<float2DReg> vec (new float2DReg(Z, X));

	if(4*n123 != sreed(data,vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

std::shared_ptr<float3DReg> VECEXT::sepRead3D(){
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
	std::shared_ptr<float3DReg> vec (new float3DReg(Z, X, Y));

	if(4*n123 != sreed("in",vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

std::shared_ptr<float3DReg> VECEXT::sepRead3D(std::string input){
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
	std::shared_ptr<float3DReg> vec (new float3DReg(Z, X, Y));

	if(4*n123 != sreed(data,vec->getVals(),4*n123))
    	seperr("Error: Cannot read data\n");

    return vec;
}

void VECEXT::copyHyper(std::shared_ptr<float2DReg> vec, const std::shared_ptr<hypercube> hyper){
    std::vector<axis> axes = hyper->getAxes();
    axes[1].n = vec->getHyper()->getAxis(2).n;
    axes[0].n = vec->getHyper()->getAxis(1).n;

    std::shared_ptr<hypercube> h (new hypercube(axes[0], axes[1]));
    vec->setHyper(h);
}

void VECEXT::add(std::shared_ptr<float2DReg> vec, const float alpha){
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    std::shared_ptr<float2D> vals = vec->_mat;    
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*vals)[ix][iz] += alpha;
        }
    }    
}

void VECEXT::set(std::shared_ptr<floatHyper> vec, const float val){
    float * pvec = vec->getVals();
    for (long long i=0; i<vec->getHyper()->getN123(); i++){
        pvec[i] = val;
    }
}

void VECEXT::set(std::shared_ptr<doubleHyper> vec, const double val){
    double * pvec = vec->getVals();
    for (long long i=0; i<vec->getHyper()->getN123(); i++){
        pvec[i] = val;
    }
}

void VECEXT::inverse(std::shared_ptr<float2DReg> vec, const float eps){
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    std::shared_ptr<float2D> vals = vec->_mat;
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*vals)[ix][iz] = 1 / ((*vals)[ix][iz] + eps);
        }
    }
}

void VECEXT::revert(std::shared_ptr<float2DReg> vec, const bool xRevert, const bool zRevert){
    unsigned int nx = vec->getHyper()->getAxis(2).n;
    unsigned int nz = vec->getHyper()->getAxis(1).n;
    std::shared_ptr<float2DReg> temp = vec->clone();

    if ((zRevert == true) && (xRevert == false)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[ix][nz-iz-1] = (*temp->_mat)[ix][iz];
            }
        }
    }

   else if ((zRevert == false) && (xRevert == true)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[nx-ix-1][iz] = (*temp->_mat)[ix][iz];
            }
        }
    }

    else if ((zRevert == true) && (xRevert == true)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[nx-ix-1][nz-iz-1] = (*temp->_mat)[ix][iz];
            }
        }
    }

    else{
        std::clog << "Nothing to revert. The input is unchanged.\n";
    }
}

void VECEXT::revert(std::shared_ptr<double2DReg> vec, const bool xRevert, const bool zRevert){
    unsigned int nx = vec->getHyper()->getAxis(2).n;
    unsigned int nz = vec->getHyper()->getAxis(1).n;
    std::shared_ptr<double2DReg> temp = vec->clone();

    if ((zRevert == true) && (xRevert == false)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[ix][nz-iz-1] = (*temp->_mat)[ix][iz];
            }
        }
    }

   else if ((zRevert == false) && (xRevert == true)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[nx-ix-1][iz] = (*temp->_mat)[ix][iz];
            }
        }
    }

    else if ((zRevert == true) && (xRevert == true)) {
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                (*vec->_mat)[nx-ix-1][nz-iz-1] = (*temp->_mat)[ix][iz];
            }
        }
    }

    else{
        std::clog << "Nothing to revert. The input is unchanged.\n";
    }
}

void VECEXT::random(std::shared_ptr<floatHyper> vec, float min, float max){
    unsigned long long n123 = vec->getHyper()->getN123();
    float * vals = vec->getVals();
    for (unsigned long long i=0; i<n123; i++){
            vals[i] += min + (max-min)*(rand() % 10000)/10000.0;
    }    
}

void VECEXT::random(std::shared_ptr<doubleHyper> vec, float min, float max){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();
    for (unsigned long long i=0; i<n123; i++){
            vals[i] += min + (max-min)*(rand() % 10000)/10000.0;
    }    
}

double VECEXT::norm2(const std::shared_ptr<floatHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    float * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i]*vals[i];
    } 
    return sum;   
}

double VECEXT::norm2(const std::shared_ptr<doubleHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i]*vals[i];
    } 
    return sum;   
}

double VECEXT::sum(const std::shared_ptr<floatHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    float * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i];
    } 
    return sum;   
}

double VECEXT::sum(const std::shared_ptr<doubleHyper> vec){
    unsigned long long n123 = vec->getHyper()->getN123();
    double * vals = vec->getVals();

    double sum = 0;
    for (unsigned long long i=0; i<n123; i++){
        sum += vals[i];
    } 
    return sum;   
}

double VECEXT::rms(const std::shared_ptr<floatHyper> vec){
    return sqrt(norm2(vec)/(vec->getHyper()->getN123()));
}

double VECEXT::rms(const std::shared_ptr<doubleHyper> vec){
    return sqrt(norm2(vec)/(vec->getHyper()->getN123()));
}

std::shared_ptr<float2DReg> VECEXT::transpose(const std::shared_ptr<float2DReg> vec){
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(X, Z));

    unsigned int nz = Z.n;
    unsigned int nx = X.n;
    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix < nx; ix++){
        for (int iz=0; iz < nz; iz++){
            (*newvals)[iz][ix] = (*vals)[ix][iz];
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::expand(const std::shared_ptr<float2DReg> vec, bool xExpand, bool zExpand) {
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nz = Z.n;
    unsigned int nx = X.n;
    std::shared_ptr<float2D> vals = vec->_mat;

    if ((zExpand == true) && (xExpand == false)) {
        Z.n = 2*nz - 1;
        Z.o *= 2;
        std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));
        std::shared_ptr<float2D> newvals = newvec->_mat;

        for (int ix=0; ix < nx; ix++){
            for (int iz=0; iz < nz; iz++){
                (*newvals)[ix][2*iz] = (*vals)[ix][iz];
            }
        }
        return newvec;
    }

    else if ((zExpand == false) && (xExpand == true)) {
        X.n = 2*nx - 1;
        X.o *= 2;
        std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));
        std::shared_ptr<float2D> newvals = newvec->_mat;

        for (int ix=0; ix < nx; ix++){
            for (int iz=0; iz < nz; iz++){
                (*newvals)[2*ix][iz] = (*vals)[ix][iz];
            }
        }
        return newvec;
    }

    else if ((zExpand == true) && (xExpand == true)) {
        X.n = 2*nx - 1;
        Z.n = 2*nz - 1;
        Z.o *= 2;
        X.o *= 2;
        std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));
        std::shared_ptr<float2D> newvals = newvec->_mat;

        for (int ix=0; ix < nx; ix++){
            for (int iz=0; iz < nz; iz++){
                (*newvals)[2*ix][2*iz] = (*vals)[ix][iz];
            }
        }
        return newvec;
    }

    else{
        std::clog << "Nothing to expand. The input is unchanged.\n";
    }
}

std::shared_ptr<float2DReg> VECEXT::interpolate(const std::shared_ptr<float2DReg> vec, bool xInterp, bool zInterp){
    
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nx = X.n;
    unsigned int nz = Z.n;
    std::shared_ptr<float2DReg> newvec;

    if ((zInterp == true) && (xInterp == false)) {
        fxTransform fx(vec->getHyper());
        std::shared_ptr<complex2DReg> cvec (new complex2DReg(fx.getRange()));
        fx.forward(false,vec,cvec);
        cvec=CVECEXT::resize(cvec,nx, nz);
        Z.d = Z.d /2;
        Z.n = 2*nz-1;
        newvec = std::make_shared<float2DReg>(Z,X);
        fx.setDomainRange(newvec->getHyper());
        fx.inverse(false,newvec,cvec);
        newvec->scale(2);
    }

    else if ((zInterp == false) && (xInterp == true)) {
        fkTransform fk(vec->getHyper());
        std::shared_ptr<complex2DReg> cvec (new complex2DReg(fk.getRange()));
        fk.forward(false,vec,cvec);
        cvec=CVECEXT::pad(cvec,floor(nx/2),floor(nx/2),0,0);
        X.d = X.d /2;
        X.n = 2*X.n-1;
        newvec = std::make_shared<float2DReg>(Z,X);
        fk.setDomainRange(newvec->getHyper());
        CVECEXT::copyHyper(cvec, fk.getRange());
        fk.inverse(false,newvec,cvec);
        newvec->scale(2);
    }

    else if ((zInterp == true) && (xInterp == true)) {
        fkTransform fk(vec->getHyper());
        std::shared_ptr<complex2DReg> cvec (new complex2DReg(fk.getRange()));
        fk.forward(false,vec,cvec);
        cvec=CVECEXT::resize(cvec,nx, nz);
        cvec=CVECEXT::pad(cvec,nx/2,nx/2,0,0);
        Z.d = Z.d /2;
        Z.n = 2*nz-1;
        X.d = X.d /2;
        X.n = 2*X.n-1;
        newvec = std::make_shared<float2DReg>(Z,X);
        fk.setDomainRange(newvec->getHyper());
        CVECEXT::copyHyper(cvec, fk.getRange());
        fk.inverse(false,newvec,cvec);
        newvec->scale(4);
    }

    else{
        newvec = vec->clone();
        std::clog << "Nothing to interpolate. The input is unchanged.\n";
    }

    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::resample(const std::shared_ptr<float2DReg> vec, bool xResamp, bool zResamp){

    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec;

    if ((zResamp == true) && (xResamp == false)) {
        Z.d = Z.d * 2;
        Z.n = ceil(Z.n / 2);
        newvec = std::make_shared<float2DReg>(Z,X);
        newvec->zero();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                (*newvec->_mat)[ix][iz] = (*vec->_mat)[ix][2*iz];
            }
        }
    }

    else if ((zResamp == false) && (xResamp == true)) {
        X.d = X.d * 2;
        X.n = ceil(X.n / 2);
        newvec = std::make_shared<float2DReg>(Z,X);
        newvec->zero();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                (*newvec->_mat)[ix][iz] = (*vec->_mat)[2*ix][iz];
            }
        }
    }

    else if ((zResamp == true) && (xResamp == true)) {
        Z.d = Z.d * 2;
        Z.n = ceil(Z.n / 2.0);
        X.d = X.d * 2;
        X.n = ceil(X.n / 2.0);
        newvec = std::make_shared<float2DReg>(Z,X);
        newvec->zero();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                (*newvec->_mat)[ix][iz] = (*vec->_mat)[2*ix][2*iz];
            }
        }
    }   

    else{
        newvec = vec->clone();
        std::clog << "Nothing to resample. The input is unchanged.\n";
    }

    return newvec;
}


std::shared_ptr<float2DReg> VECEXT::resize(const std::shared_ptr<float2DReg> vec, unsigned int nx, unsigned int nz){
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    X.n = nx;
    Z.n = nz;
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix < std::min(nx,VECEXT::getnx(vec)); ix++){
        for (int iz=0; iz < std::min(nz,VECEXT::getnz(vec)); iz++){
            (*newvals)[ix][iz] = (*vals)[ix][iz];
        }
    }
    return newvec;
}

std::shared_ptr<float1DReg> VECEXT::pad(const std::shared_ptr<float1DReg> vec, unsigned int nz_start, unsigned int nz_end){

    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nz = Z.n;
    Z.n += nz_start + nz_end;
    Z.o = Z.o - nz_start*Z.d;

    std::shared_ptr<float1DReg> newvec (new float1DReg(Z));
    newvec->zero();

    std::shared_ptr<float1D> vals = vec->_mat;
    std::shared_ptr<float1D> newvals = newvec->_mat;

    for (int iz=nz_start; iz<nz_start+nz; iz++){
            (*newvals)[iz] = (*vals)[iz-nz_start];
    }
    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::pad(const std::shared_ptr<float2DReg> vec, unsigned int nx_start, unsigned int nx_end,
                                unsigned int nz_start, unsigned int nz_end){
    
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nx = X.n;
    unsigned int nz = Z.n;
    X.n += nx_start + nx_end;
    Z.n += nz_start + nz_end;
    X.o = X.o - nx_start*X.d;
    Z.o = Z.o - nz_start*Z.d;

    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));
    newvec->zero();

    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=nx_start; ix<nx_start+nx; ix++){
        for (int iz=nz_start; iz<nz_start+nz; iz++){
            (*newvals)[ix][iz] = (*vals)[ix-nx_start][iz-nz_start];
        }
    }
    return newvec;
}

std::shared_ptr<float1DReg> VECEXT::cut(const std::shared_ptr<float1DReg> vec, unsigned int nz_start, unsigned int nz_end){
    
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nz = Z.n;
    Z.n = Z.n - nz_start - nz_end;
    Z.o = Z.o + nz_start*Z.d;

    std::shared_ptr<float1DReg> newvec (new float1DReg(Z));
    newvec->zero();

    std::shared_ptr<float1D> vals = vec->_mat;
    std::shared_ptr<float1D> newvals = newvec->_mat;

    for (int iz=0; iz<Z.n; iz++){
        (*newvals)[iz] = (*vals)[iz+nz_start];
    }
    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::cut(const std::shared_ptr<float2DReg> vec, unsigned int nx_start, unsigned int nx_end,
                                unsigned int nz_start, unsigned int nz_end){
    
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nx = X.n;
    unsigned int nz = Z.n;
    X.n = X.n - nx_start - nx_end;
    Z.n = Z.n - nz_start - nz_end;
    X.o = X.o + nx_start*X.d;
    Z.o = Z.o + nz_start*Z.d;

    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));
    newvec->zero();

    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*newvals)[ix][iz] = (*vals)[ix+nx_start][iz+nz_start];
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::reshape(const std::shared_ptr<float2DReg> vec1, const std::shared_ptr<float2DReg> vec2){

    axis X1 = vec1->getHyper()->getAxis(2);
    axis Z1 = vec1->getHyper()->getAxis(1);
    axis X2 = vec2->getHyper()->getAxis(2);
    axis Z2 = vec2->getHyper()->getAxis(1);

    if((X1.d != X2.d) || (Z1.d != Z2.d))
        throw std::logic_error("Sampling mismatch, vector cannot be reshaped\n");
    
    if( (abs(round((X1.o - X2.o)/X1.d) - (X1.o - X2.o)/X1.d) > 1e-05) || (abs(round((Z1.o - Z2.o)/Z1.d) - (Z1.o - Z2.o)/Z1.d) > 1e-05) )
        throw std::logic_error("Origins have wrong multiplicity, vector cannot be reshaped\n");
    
    std::shared_ptr<float2DReg> newvec = vec2->clone();

    int min2_nx, min2_nz, max2_nx, max2_nz, min1_nx, max1_nx, min1_nz, max1_nz;
    min1_nx = (X1.o-X2.o)/X1.d;
    min2_nx = std::max(0, min1_nx);
    max1_nx = X2.n - (X2.o+X2.n*X2.d - X1.o-X1.n*X1.d)/X1.d;
    max2_nx = std::min(X2.n, max1_nx);
    min1_nz = (Z1.o-Z2.o)/Z1.d;
    min2_nz = std::max(0, min1_nz);
    max1_nz = Z2.n - (Z2.o+Z2.n*Z2.d - Z1.o-Z1.n*Z1.d)/Z1.d;
    max2_nz = std::min(Z2.n, max1_nz);

    for (int ix=min2_nx; ix<max2_nx; ix++){
        for (int iz=min2_nz; iz<max2_nz; iz++){
            (*newvec->_mat)[ix][iz] = (*vec1->_mat)[ix-min1_nx][iz-min1_nz];
        }
    }
    return newvec;
}

std::shared_ptr<float1DReg> VECEXT::stackx(const std::shared_ptr<float2DReg> vec){
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    axis Z = vec->getHyper()->getAxis(1);
    std::shared_ptr<float1DReg> newvec (new float1DReg(Z));

    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float1D> newvals = newvec->_mat;

    for (int ix=0; ix < nx; ix++){
        for (int iz=0; iz < nz; iz++){
            (*newvals)[iz] += (*vals)[ix][iz];
        }
    }
    newvec->scale(1.0/nx);

    return newvec;
}

std::shared_ptr<float1DReg> VECEXT::stackz(const std::shared_ptr<float2DReg> vec){
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    axis Z = vec->getHyper()->getAxis(1);
    std::shared_ptr<float1DReg> newvec (new float1DReg(Z));

    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<float1D> newvals = newvec->_mat;

    for (int ix=0; ix < nx; ix++){
        for (int iz=0; iz < nz; iz++){
            (*newvals)[ix] += (*vals)[ix][iz];
        }
    }
    newvec->scale(1.0/nz);

    return newvec;
}

void VECEXT::lpSinc(std::shared_ptr<float2DReg> vec, float wc, int N, float a){

    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    if (2*N+1>nz)
        throw std::logic_error("filter longer than data\n");

    std::shared_ptr<float2D> vals = vec->_mat;

    for (int ix=0; ix < nx; ix++){

        // store filter coefficients and needed portion of the data in temporary vectors
        std::vector <float> temp, filt;
        float num;
        
        // apply filtering in the roll-on part
        for (int i=0; i < N; i++){
            temp.push_back((*vals)[ix][i]);
            (*vals)[ix][i]=0;
            num = wc*sinc(wc*M_PI*(i-N))*(a +(1-a)*cos(M_PI*(i-N)/N));
            filt.push_back(num);
        }

        for (int i=N; i < 2*N+1; i++){			
            temp.push_back((*vals)[ix][i]);
            (*vals)[ix][i]=0;
            num = wc*sinc(wc*M_PI*(i-N))*(a +(1-a)*cos(M_PI*(i-N)/N));
            filt.push_back(num);
            for (int j=0; j<=i; j++){
                (*vals)[ix][i-N] = (*vals)[ix][i-N] + temp[i-j]*filt[j];
            }
        }

        // apply filtering to the rest of the data
        for (int i=2*N+1; i<nz; i++){
            temp.erase(temp.begin());
            temp.push_back((*vals)[ix][i]);
            (*vals)[ix][i-N]=0;
            for (int j=0; j<2*N+1; j++){
                (*vals)[ix][i-N] = (*vals)[ix][i-N] + temp[2*N-j]*filt[j];
            }
        }

        // apply filtering in the roll-off part
        for (int i=nz; i < N+nz; i++){
            temp.erase(temp.begin());
            filt.erase(filt.begin());
            (*vals)[ix][i-N]=0;
            for (int j=0; j<temp.size(); j++){
                (*vals)[ix][i-N] = (*vals)[ix][i-N] + temp[temp.size()-j]*filt[j];
            } 
        }
    }
}

void VECEXT::iirButterworth(std::shared_ptr<float1DReg> vec, float wc, bool lp, unsigned int ho){

    if (lp == false) wc = 1-wc;

    // Initialization
    unsigned int nz = VECEXT::getnz(vec);
    std::shared_ptr<float1D> vals = vec->_mat;

    int M = 2*ho;
    float Omega_c = tan(M_PI*wc/2);
    float a0, a1, a2, b1, b2, c;
    float xi, xi1, xi2;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        c=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_c + Omega_c*Omega_c;
        a0=Omega_c*Omega_c/c;
        a1 = 2 * a0;
        a2 = a0;
        b1=2*(Omega_c*Omega_c-1)/c;
        b2=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_c +Omega_c*Omega_c)/c;

        // deduce filter coefficients for High-Pass
        if (lp == false){
            a0 = a0;
            a1 = - a1;
            a2 = a2;
            b1 = - b1;
            b2 = b2;
        }

        // apply filter recursively
        xi = (*vals)[1];
        xi1 = (*vals)[0];
        xi2 = 0;
        (*vals)[0] = a0*xi1;
        (*vals)[1] = a0*xi + a1*xi1 - b1*(*vals)[0];

        for (int iz=2; iz<nz; iz++){
            xi2 = xi1;
            xi1 = xi;			
            xi = (*vals)[iz];
            (*vals)[iz] = a0*xi + a1*xi1 + a2*xi2 
                    - b1*(*vals)[iz-1] -b2*(*vals)[iz-2];
        }
    }
}

void VECEXT::iirButterworth(std::shared_ptr<float2DReg> vec, float wc, bool lp, unsigned int ho){

    if (lp == false) wc = 1-wc;

    // Initialization
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    std::shared_ptr<float2D> vals = vec->_mat;

    int M = 2*ho;
    float Omega_c = tan(M_PI*wc/2);
    float a0, a1, a2, b1, b2, c;
    float xi, xi1, xi2;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        c=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_c + Omega_c*Omega_c;
        a0=Omega_c*Omega_c/c;
        a1 = 2 * a0;
        a2 = a0;
        b1=2*(Omega_c*Omega_c-1)/c;
        b2=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_c +Omega_c*Omega_c)/c;

        // deduce filter coefficients for High-Pass
        if (lp == false){
            a0 = a0;
            a1 = - a1;
            a2 = a2;
            b1 = - b1;
            b2 = b2;
        }

        // apply filter recursively
        for (int ix=0; ix < nx; ix++){
            xi = (*vals)[ix][1];
            xi1 = (*vals)[ix][0];
            xi2 = 0;
            (*vals)[ix][0] = a0*xi1;
            (*vals)[ix][1] = a0*xi + a1*xi1 - b1*(*vals)[ix][0];

            for (int iz=2; iz<nz; iz++){
                xi2 = xi1;
                xi1 = xi;			
                xi = (*vals)[ix][iz];
                (*vals)[ix][iz] = a0*xi + a1*xi1 + a2*xi2 
                        - b1*(*vals)[ix][iz-1] -b2*(*vals)[ix][iz-2];
            }
        }
    }
}

void VECEXT::iirBpButterworth(std::shared_ptr<float1DReg> vec, float wcl, float wch, unsigned int ho){
    
    // Initialization
    unsigned int nz = VECEXT::getnz(vec);
    std::shared_ptr<float1D> vals = vec->_mat;

    int M = 2*ho;
    float Omega_ch = tan(M_PI*wch/2);
    float Omega_cl = tan(M_PI*(1-wcl)/2);
    float a0h, a1h, a2h, b1h, b2h, ch;
    float a0l, a1l, a2l, b1l, b2l, cl;
    float a0, a1, a2, a3, a4, b1, b2, b3, b4;
    float xi, xi1, xi2, xi3, xi4;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        ch=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch + Omega_ch*Omega_ch;
        a0h=Omega_ch*Omega_ch/ch;
        a1h = 2 * a0h;
        a2h = a0h;
        b1h=2*(Omega_ch*Omega_ch-1)/ch;
        b2h=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch +Omega_ch*Omega_ch)/ch;

        // compute filter coefficients for High-Pass
        cl=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl + Omega_cl*Omega_cl;
        a0l=Omega_cl*Omega_cl/cl;
        a1l = - 2 * a0l;
        a2l = a0l;
        b1l=-2*(Omega_cl*Omega_cl-1)/cl;
        b2l=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl +Omega_cl*Omega_cl)/cl;

        // deduce filter coefficients for Band-Pass
        a0 = a0h*a0l;
        a1 = a0h*a1l + a0l*a1h;
        a2 = a0h*a2l + a1h*a1l + a2h*a0l;
        a3 = a1h*a2l + a2h*a1l;
        a4 = a2h*a2l;
        b1 = b1l + b1h;
        b2 = b2l + b1h*b1l + b2h;
        b3 = b1h*b2l + b2h*b1l;
        b4 = b2h*b2l;
        
        // apply filter recursively
        xi = (*vals)[3];
        xi1 = (*vals)[2];
        xi2 = (*vals)[1];
        xi3 = (*vals)[0];
        xi4 = 0;
        (*vals)[0] = a0*xi3;
        (*vals)[1] = a0*xi2 + a1*xi3 - b1*(*vals)[0];
        (*vals)[2] = a0*xi1 + a1*xi2 + a2*xi3 - b1*(*vals)[1] - b2*(*vals)[0];
        (*vals)[3] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 - b1*(*vals)[2] - b2*(*vals)[1] - b3*(*vals)[0];

        for (int iz=4; iz<nz; iz++){
            xi4 = xi3;
            xi3 = xi2;
            xi2 = xi1;
            xi1 = xi;			
            xi = (*vals)[iz];
            (*vals)[iz] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 + a4*xi4
                    - b1*(*vals)[iz-1] -b2*(*vals)[iz-2] -b3*(*vals)[iz-3] -b4*(*vals)[iz-4];
        }
    }
}

void VECEXT::iirBpButterworth(std::shared_ptr<float2DReg> vec, float wcl, float wch, unsigned int ho){
    
    // Initialization
    unsigned int nz = VECEXT::getnz(vec);
    unsigned int nx = VECEXT::getnx(vec);
    std::shared_ptr<float2D> vals = vec->_mat;

    int M = 2*ho;
    float Omega_ch = tan(M_PI*wch/2);
    float Omega_cl = tan(M_PI*(1-wcl)/2);
    float a0h, a1h, a2h, b1h, b2h, ch;
    float a0l, a1l, a2l, b1l, b2l, cl;
    float a0, a1, a2, a3, a4, b1, b2, b3, b4;
    float xi, xi1, xi2, xi3, xi4;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        ch=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch + Omega_ch*Omega_ch;
        a0h=Omega_ch*Omega_ch/ch;
        a1h = 2 * a0h;
        a2h = a0h;
        b1h=2*(Omega_ch*Omega_ch-1)/ch;
        b2h=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch +Omega_ch*Omega_ch)/ch;

        // compute filter coefficients for High-Pass
        cl=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl + Omega_cl*Omega_cl;
        a0l=Omega_cl*Omega_cl/cl;
        a1l = - 2 * a0l;
        a2l = a0l;
        b1l=-2*(Omega_cl*Omega_cl-1)/cl;
        b2l=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl +Omega_cl*Omega_cl)/cl;

        // deduce filter coefficients for Band-Pass
        a0 = a0h*a0l;
        a1 = a0h*a1l + a0l*a1h;
        a2 = a0h*a2l + a1h*a1l + a2h*a0l;
        a3 = a1h*a2l + a2h*a1l;
        a4 = a2h*a2l;
        b1 = b1l + b1h;
        b2 = b2l + b1h*b1l + b2h;
        b3 = b1h*b2l + b2h*b1l;
        b4 = b2h*b2l;
        
        // apply filter recursively
        for (int ix=0; ix < nx; ix++){
            xi = (*vals)[ix][3];
            xi1 = (*vals)[ix][2];
            xi2 = (*vals)[ix][1];
            xi3 = (*vals)[ix][0];
            xi4 = 0;
            (*vals)[ix][0] = a0*xi3;
            (*vals)[ix][1] = a0*xi2 + a1*xi3 - b1*(*vals)[ix][0];
            (*vals)[ix][2] = a0*xi1 + a1*xi2 + a2*xi3 - b1*(*vals)[ix][1] - b2*(*vals)[ix][0];
            (*vals)[ix][3] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 - b1*(*vals)[ix][2] - b2*(*vals)[ix][1] - b3*(*vals)[ix][0];

            for (int iz=4; iz<nz; iz++){
                xi4 = xi3;
                xi3 = xi2;
                xi2 = xi1;
                xi1 = xi;			
                xi = (*vals)[ix][iz];
                (*vals)[ix][iz] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 + a4*xi4
                        - b1*(*vals)[ix][iz-1] -b2*(*vals)[ix][iz-2] -b3*(*vals)[ix][iz-3] -b4*(*vals)[ix][iz-4];
            }
        }
    }
}

void VECEXT::fkFilter(std::shared_ptr<float2DReg> vec, float kLow, float kHigh, float taper){
    if ((abs(kLow)>0.5) || (abs(kHigh)>0.5) || (kLow>kHigh))
        throw std::logic_error("The cutoff wavenumbers must obey -0.5 < kLow < kHigh < 0.5");

    fkTransform fk(vec->getHyper());
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    Z.n = Z.n /2 +1;
    X.d = fk.dk();
    X.o = fk.ok();
    Z.d = fk.df();
    Z.o = 0;
    std::shared_ptr<complex2DReg> cvec (new complex2DReg(Z, X));
    fk.forward(false,vec,cvec);

    int nx1 = (0.5+kLow-taper)*X.n;
    nx1 = std::max(nx1, 0);
    int nx2 = (0.5+kLow)*X.n;
    int nx3 = std::ceil((0.5+kHigh)*X.n);
    int nx4 = (0.5+kHigh+taper)*X.n;
    nx4 = std::min(nx4, X.n);
    
    for (int ix=0; ix<nx1; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*cvec->_mat)[ix][iz] = 0;
        }
    }
    
    for (int ix=nx1; ix<nx2; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*cvec->_mat)[ix][iz] *= (1 + std::cos(M_PI*(kLow + 0.5 - (1.0*ix)/X.n)/taper))/2.0;
        }
    }
    
    for (int ix=nx3; ix<nx4; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*cvec->_mat)[ix][iz] *= (1 + std::cos(M_PI*((1.0*ix)/X.n - 0.5 - kHigh)/taper))/2.0;
        }
    }
    
    for (int ix=nx4; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*cvec->_mat)[ix][iz] = 0;
        }
    }
    
    fk.inverse(false,vec,cvec);
}

std::shared_ptr<float1DReg> VECEXT::zero_phase(const std::shared_ptr<float1DReg> vec){
    
    axis Z = vec->getHyper()->getAxis(1);
    axis X = Z;
    X.n = 1;
    
    // Make a new vector with odd length (fake 2D to be able to use fxtranform operator)
    std::shared_ptr<float2DReg> vec0;
	if (2*(Z.n/2)==Z.n) {
		Z.n+=1;
        vec0 = std::make_shared<float2DReg>(Z, X);

        for (int iz=0; iz<Z.n-1; iz++){
        vec0->getVals()[iz] = vec->getVals()[iz];
        }
        vec0->getVals()[Z.n-1] = 0;
    }
	
    else {
        vec0 = std::make_shared<float2DReg>(Z, X);

        for (int iz=0; iz<Z.n; iz++){
        vec0->getVals()[iz] = vec->getVals()[iz];
        }
    }

    // Fourier transform the vector
    fxTransform fx(vec0->getHyper());
    axis Z2;
    Z2.d = fx.df();
    Z2.o = 0;
    Z2.n = Z.n /2 + 1;
	std::shared_ptr<complex2DReg> cvec0 (new complex2DReg(Z2, X));
	fx.forward(false,vec0,cvec0);

    int N = Z2.n;
    double dw = 2*M_PI/Z.n;
    double r;
    std::complex<double> z (0,0);

    // Modify the FT by setting a linear phase while keeping the same power spectrum
	for (int i=0; i<N; i++){

        r = sqrt(norm(cvec0->getVals()[i]));

        // add linear phase shift (* exp (-j.omega.dt.(Nt-1)/2)))
        z.real(r*cos(-i*dw*(Z.n - 1)/2));
        z.imag(r*sin(-i*dw*(Z.n - 1)/2));

        cvec0->getVals()[i] = z;
	}

    // inverse FT and modify the original vector
    fx.inverse(false,vec0,cvec0);
    Z.o = (int)(-Z.n/2) * Z.d; 
    std::shared_ptr<float1DReg> vec1 (new float1DReg(Z));
    for (int iz=0; iz<Z.n; iz++){
        vec1->getVals()[iz] = vec0->getVals()[iz];
    }
    
    return vec1;
}

std::shared_ptr<float1DReg> VECEXT::minimum_phase(const std::shared_ptr<float1DReg> vec, float eps){
    
    axis Z = vec->getHyper()->getAxis(1);
    axis X = Z;
    X.n = 1;
    
    // Make a new vector with odd length (fake 2D to be able to use fxtranform operator)
    std::shared_ptr<float2DReg> vec0;
	if (2*(Z.n/2)==Z.n) {
		Z.n+=1;
        vec0 = std::make_shared<float2DReg>(Z, X);

        for (int iz=0; iz<Z.n-1; iz++){
        vec0->getVals()[iz] = vec->getVals()[iz];
        }
        vec0->getVals()[Z.n-1] = 0;
    }
	
    else {
        vec0 = std::make_shared<float2DReg>(Z, X);

        for (int iz=0; iz<Z.n; iz++){
        vec0->getVals()[iz] = vec->getVals()[iz];
        }
    }

    // Fourier transform the vector
    fxTransform fx(vec0->getHyper());
    axis Z2;
    Z2.d = fx.df();
    Z2.o = 0;
    Z2.n = Z.n /2 + 1;
	std::shared_ptr<complex2DReg> cvec0 (new complex2DReg(Z2, X));
	fx.forward(false,vec0,cvec0);

    int N = Z2.n;
    double r;
    std::complex<double> z (0,0);

    // Modify the FT by making it real and equal to the logarithm of the power spectrum
	for (int i=0; i<N; i++){

        r = sqrt(norm(cvec0->getVals()[i]));

        // set FT = log(r+eps)
        z.real(log(r+eps));
        z.imag(0);

        cvec0->getVals()[i] = z;
	}

    // inverse FT and keep the causal part
    fx.inverse(false,vec0,cvec0);
    for (int iz=1; iz<Z.n/2+1; iz++){
        vec0->getVals()[iz] *= 2;
    }
    for (int iz=Z.n/2+1; iz<Z.n; iz++){
        vec0->getVals()[iz] = 0;
    }

    // FT again and take the exponent
    fx.forward(false,vec0,cvec0);
    for (int iz=0; iz<N; iz++){
        cvec0->getVals()[iz] = exp(cvec0->getVals()[iz]);
    }

    // inverse FT
    fx.inverse(false,vec0,cvec0);
    Z.o = 0;
    std::shared_ptr<float1DReg> vec1 (new float1DReg(Z));
    for (int iz=0; iz<Z.n; iz++){
        vec1->getVals()[iz] = vec0->getVals()[iz];
    }
    
    return vec1;
}

std::shared_ptr<float1DReg> VECEXT::alphaTrim(const std::shared_ptr<float1DReg> vec, unsigned int halfLength, float alpha){
    if ((alpha>1) || (alpha<0))
        throw std::logic_error("alpha must obey 0 <= alpha <= 1");

    if (vec->getHyper()->getN123()<2*halfLength+1)
        throw std::logic_error("Input vector too short of filter length too long");

    std::shared_ptr<float1DReg> newvec = vec->clone();
    newvec->zero();

    int N = vec->getHyper()->getN123();
    int trim = floor(alpha*halfLength);
    float avg=0;

    // Extend the input vector symmetrically wrt to the boundaries
    std::shared_ptr<float1DReg> vec_ex = pad(vec, halfLength, halfLength);

    for (int i=0; i<halfLength; i++){
        vec_ex->getVals()[i] = vec->getVals()[halfLength-i-1];
        vec_ex->getVals()[N+halfLength+i] = vec->getVals()[N-i-1];
    }

    std::vector<float> kernel(2*halfLength+1);

    // Apply the smoothing
    for (int i=0; i<N; i++){

        // Copy values to the kernel
        for (int j=0; j<2*halfLength+1; j++){
            kernel[j] = vec_ex->getVals()[i+j];
        }

        // Sort the values in the kernel
        std::sort(kernel.begin(), kernel.end());

        // Compute the average after dropping the trimmed values
        avg = 0;
        for (int j=trim; j<2*halfLength+1-trim; j++){
            avg += kernel[j];
        }
        avg = avg / (2*halfLength+1-2*trim);

        // Copy the value to the output
        newvec->getVals()[i] = avg;

    }

    return newvec;
}

std::shared_ptr<float2DReg> VECEXT::alphaTrim(const std::shared_ptr<float2DReg> vec, unsigned int halfLengthx, unsigned int halfLengthz, float alpha){
    if ((alpha>1) || (alpha<0))
        throw std::logic_error("alpha must obey 0 <= alpha <= 1");

    if (vec->getHyper()->getAxis(1).n<2*halfLengthz+1)
        throw std::logic_error("Input vector too short of filter length too long in the Z direction");

    if (vec->getHyper()->getAxis(2).n<2*halfLengthx+1)
        throw std::logic_error("Input vector too short of filter length too long in the X direction");

    std::shared_ptr<float2DReg> newvec = vec->clone();
    newvec->zero();

    int nz = vec->getHyper()->getAxis(1).n;
    int nx = vec->getHyper()->getAxis(2).n;
    int trim = floor(alpha*((2*halfLengthz+1)*(2*halfLengthx+1)-1)/2 );
    float avg=0;

    // Extend the input vector symmetrically wrt to the boundaries
    std::shared_ptr<float2DReg> vec_ex = pad(vec, halfLengthx, halfLengthx, halfLengthz, halfLengthz);

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<halfLengthz; iz++){
            (*vec_ex->_mat)[ix+halfLengthx][iz] = (*vec->_mat)[ix][halfLengthz-iz-1];
            (*vec_ex->_mat)[ix+halfLengthx][nz+halfLengthz+iz] = (*vec->_mat)[ix][nz-iz-1];
        }
    }

    for (int iz=0; iz<nz+2*halfLengthz; iz++){
        for (int ix=0; ix<halfLengthx; ix++){
            (*vec_ex->_mat)[ix][iz] = (*vec_ex->_mat)[2*halfLengthx-ix-1][iz];
            (*vec_ex->_mat)[nx+halfLengthx+ix][iz] = (*vec_ex->_mat)[nx+halfLengthx-ix-1][iz];
        }
    }

    std::vector<float> kernel((2*halfLengthx+1)*(2*halfLengthz+1));

    // Apply the smoothing
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz; iz++){

            // Copy values to the kernel
            for (int jx=0; jx<2*halfLengthx+1; jx++){
                for (int jz=0; jz<2*halfLengthz+1; jz++){
                    kernel[jx*(2*halfLengthz+1)+jz]=(*vec_ex->_mat)[ix+jx][iz+jz];
                }
            }

            // Sort the values in the kernel
            std::sort(kernel.begin(), kernel.end());

            // Compute the average after dropping the trimmed values
            avg = 0;
            for (int j=trim; j<(2*halfLengthx+1)*(2*halfLengthz+1)-trim; j++){
                avg += kernel[j];
            }
            avg = avg / ((2*halfLengthx+1)*(2*halfLengthz+1)-2*trim);

            // Copy the value to the output
            (*newvec->_mat)[ix][iz] = avg;
        }
    }

    return newvec;
}
