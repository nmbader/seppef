#include <time.h>
#include <chrono>

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "floatHyperExt.h"
#include "float2DReg.h"
#include "functions.h"
#include "zeroPhaseFiltering2D.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program filter \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// list of available filters
std::vector<std::string> list = {"butterworth", "sinc"};

// read in the filter type and other parameters
    std::string type;
    bool zero_phase;
    float lwc, hwc, alpha;
    int hl, ho;
    MB::readParam(argc, argv, "type", type, "sinc");
    MB::readParam(argc, argv, "zero_phase", zero_phase, false);
    MB::readParam(argc, argv, "low_cutoff", lwc, 0);
    MB::readParam(argc, argv, "high_cutoff", hwc, 1);
    MB::readParam(argc, argv, "half_length", hl, 20);
    MB::readParam(argc, argv, "half_order", ho, 2);
    MB::readParam(argc, argv, "alpha", alpha, 0.5);

    std::shared_ptr<float2DReg> out (new float2DReg(dat->getHyper()));

    axis Z, Zf;
    Z = dat->getHyper()->getAxis(1);

// construct and apply the filters
    if (type=="sinc"){
        std::shared_ptr<float1DReg> lpfilter = MB::sincWavelet(hwc,hl,alpha);
        std::shared_ptr<float1DReg> lpfilter2 = MB::sincWavelet(lwc,hl,alpha);
        Zf = lpfilter->getHyper()->getAxis(1);
        Zf.d = Z.d;
        Zf.o = Z.o - Zf.d*hl;
        std::shared_ptr<hypercube> hyperf (new hypercube(Zf));
        lpfilter->setHyper(hyperf);
        lpfilter2->setHyper(hyperf);
        std::shared_ptr<float1DReg> hpfilter = lpfilter2->clone();
        hpfilter->scale(-1);
        (*hpfilter->_mat)[hl] += 1;

        zeroPhaseFiltering2D Fl(lpfilter,dat->getHyper(),dat->getHyper());
        zeroPhaseFiltering2D Fh(hpfilter,dat->getHyper(),dat->getHyper());

        if (lwc > 0){
            Fh.forward(false,dat,out);
            dat = out->clone();
        }

        if (hwc<1){
            Fl.forward(false,dat,out);
        }
    }

    else if ((type=="butterworth") && (zero_phase == false)){
        if (lwc == 0) {
            VECEXT::iirButterworth(dat,hwc,true,ho);
        }
        else if (hwc == 1){
            VECEXT::iirButterworth(dat,lwc,false,ho);
        }
        else{
            VECEXT::iirBpButterworth(dat,lwc,hwc,ho);
        }
        out = dat;
    }

    else if ((type=="butterworth") && (zero_phase == true)){

        std::shared_ptr<float1DReg> filter (new float1DReg(Z));
        filter->getVals()[0] = 1;

        if (lwc == 0) {
            VECEXT::iirButterworth(filter,hwc,true,ho);
        }
        else if (hwc == 1){
            VECEXT::iirButterworth(filter,lwc,false,ho);
        }
        else{
            VECEXT::iirBpButterworth(filter,lwc,hwc,ho);
        }
        filter = VECEXT::zero_phase(filter);
        Zf.n = filter->getHyper()->getAxis(1).n;
        
        int cut = Zf.n -2*hl -1;
        cut = std::max(0, cut/2);
        filter = VECEXT::cut(filter, cut, cut);

        Zf.n = filter->getHyper()->getAxis(1).n;
        Zf.d = Z.d;
        Zf.o = Z.o - Zf.d*hl;
        std::shared_ptr<hypercube> hyperf (new hypercube(Zf));
        filter->setHyper(hyperf);

        zeroPhaseFiltering2D Fil(filter,dat->getHyper(),dat->getHyper());
        Fil.forward(false,dat,out);
    }

    else{
        throw std::logic_error("The requested filter is not implemented.\n");
    }


// write the filtered data to disk
    VECEXT::sepWrite(out, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program filter; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}