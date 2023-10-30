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
    
    std::clog << "Starting program diverge to apply spherical divergence (t/tref)^power \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read dimensions form history file
    int n2, n3;
    if(0 == hetch((char*)"n3",(char*)"d",&n3))
    	n3=0;
    if(0 == hetch((char*)"n2",(char*)"d",&n2))
    	n2=0;

// read in the input data base on the dimension above
    std::shared_ptr<floatHyper> dat;
    if (n2==0)
        dat = VECEXT::sepRead1D(); // 1D vector
    else if (n3==0)
        dat = VECEXT::sepRead(); // 2D vector

    else
        dat = VECEXT::sepRead3D(); // 3D vector
    

// read in parameters
    float tref, tpower;
    bool reverse;

    if(0 == getch("tref",(char*)"f",&tref))
        tref=0.1;

    if(0 == getch("tpower",(char*)"f",&tpower))
        tpower=2.0;

    MB::readParam(argc, argv, "reverse", reverse, false);

// apply (or back-off) spherical divergence to data
    axis Z = dat->getHyper()->getAxis(1);
    long long n123 = dat->getHyper()->getN123();
    float * ptr = dat->getVals();
    float time;
    if (reverse == false){
        for (long long i=0; i<n123; i++){
            time = (i - floor(i/Z.n)*Z.n)*Z.d-Z.o;
            ptr[i] *= pow((time/tref),tpower);
        }
    }
    else{
        for (long long i=0; i<n123; i++){
            time = (i - floor(i/Z.n)*Z.n)*Z.d-Z.o;
            ptr[i] *= pow((time/tref),-tpower);
        }
    }

// write the output to disk
    VECEXT::sepWrite(dat, output);
	
    return 0;
}