#include "float2DReg.h"
#include "floatHyperExt.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program resamp to resample data \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in parameters
    bool xresamp, zresamp;

    MB::readParam(argc, argv, "xresamp", xresamp, true);
    MB::readParam(argc, argv, "zresamp", zresamp, true);

// resample data
    dat = VECEXT::resample(dat, xresamp, zresamp);

// write the output to disk
    VECEXT::sepWrite(dat, output);
	
    return 0;
}