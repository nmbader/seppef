#include "float1DReg.h"
#include "float2DReg.h"
#include "float3DReg.h"
#include "floatHyperExt.h"
#include "axis.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program duplicate \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read dimensions form history file
    int n2, n3;
    if(0 == hetch((char*)"n3",(char*)"d",&n3))
    	n3=0;
    if(0 == hetch((char*)"n2",(char*)"d",&n2))
    	n2=0;

// read in parameters
    int nb_copies;
    if(0 == getch("nb_copies",(char*)"d",&nb_copies))
        nb_copies=1;

// read in the input data base on the dimension above
    std::shared_ptr<floatHyper> dat;
    std::shared_ptr<floatHyper> out;
    axis Z, X, T;
    if (n2==0){
        dat = VECEXT::sepRead1D(); // 1D vector
        Z = dat->getHyper()->getAxis(1);
        X.n = nb_copies;
        out = std::make_shared<float2DReg> (Z,X);
    }
    else if (n3==0){
        dat = VECEXT::sepRead(); // 2D vector
        Z = dat->getHyper()->getAxis(1);
        X = dat->getHyper()->getAxis(2);
        T.n = nb_copies;
        out = std::make_shared<float3DReg> (Z,X,T);
    }
    else {
        dat = VECEXT::sepRead3D(); // 3D vector
        Z = dat->getHyper()->getAxis(1);
        X = dat->getHyper()->getAxis(2);
        T = dat->getHyper()->getAxis(3);
        T.n = T.n * nb_copies;
        out = std::make_shared<float3DReg> (Z,X,T);
    }
    

// create copies of the data
    long long n123 = dat->getHyper()->getN123();
    float * ptr = dat->getVals();
    float * ptr2 = out->getVals();
    for (int i=0; i<nb_copies; i++){
        for (long long j=0; j<n123; j++){
            ptr2[n123*i+j] = ptr[j];
        }
    }


// write the output to disk
    VECEXT::sepWrite(out, output);
	
    return 0;
}