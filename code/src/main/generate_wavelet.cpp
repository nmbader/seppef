#include <time.h>
#include <chrono>

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "floatHyperExt.h"
#include "float1DReg.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in parameters
    float wc, dz, amp, lwc, hwc, alpha, sigma, epsilon, rho;
    int nz, shift, ho;
    std::string type, phase;
    
    if(0 == getch("wc",(char*)"f",&wc)) // central freq for ricker
        wc=0.1;

    if(0 == getch("sigma",(char*)"f",&sigma)) // std dev in nb of samples for Gaussian
        sigma=10;

    if(0 == getch("low_cutoff",(char*)"f",&lwc)) // cutoff for sinc and butterworth
        lwc=0.0;

    if(0 == getch("high_cutoff",(char*)"f",&hwc))
        hwc=1.0;

    if(0 == getch("alpha",(char*)"f",&alpha)) // used for cosine window for sinc
        alpha=0.5;

    if(0 == getch("half_order",(char*)"d",&ho)) // half order used for butterworth
        ho=2;

    if(0 == getch("epsilon",(char*)"f",&epsilon)) // used in Kolgomoroff factorization for minimum phase
        epsilon=1e-07;

    if(0 == getch("rho",(char*)"f",&rho)) // used in leaky integration
        rho=1e-07;

    if(0 == getch("dz",(char*)"f",&dz))
        dz=0.004;    

    if(0 == getch("amp",(char*)"f",&amp))
        amp=1;

    if(0 == getch("nz",(char*)"d",&nz))
        nz=100;

    if(0 == getch("shift",(char*)"d",&shift))
        shift=0;

    MB::readParam(argc, argv, "type", type, "ricker");
    MB::readParam(argc, argv, "phase", phase, "default");


// construct the wavelet
    std::shared_ptr<float1DReg> dat;
    int half_length = floor(nz/2);
    int length = 2*half_length+1;

    if (type=="ricker")
	    dat = MB::rickerWavelet(wc, half_length);
    else if (type=="gaussian")
        dat=MB::gaussianWavelet(sigma, half_length);
    else if (type=="sinc"){
        std::shared_ptr<float1DReg> lp = MB::sincWavelet(hwc, half_length, alpha);
        std::shared_ptr<float1DReg> hp = MB::sincWavelet(lwc, half_length, alpha);
        hp->scale(-1);
        hp->getVals()[half_length] += 1;
        if (lwc==0)
            dat=lp->clone();
        else if (hwc==1)
            dat=hp->clone();
        else{
            dat=lp->clone();
            dat->zero();
            for (int i=half_length; i<length; i++){
                for (int j=0; j<=i; j++){
                    dat->getVals()[i-half_length] += lp->getVals()[i-j]*hp->getVals()[j];
                }
                dat->getVals()[length-1+half_length-i] = dat->getVals()[i-half_length];
            }
        }
    }
    else if (type=="butterworth"){
        dat=std::make_shared<float1DReg> (length);
        dat->zero();
        dat->getVals()[0] = 1;
        if (lwc==0)
            VECEXT::iirButterworth(dat,hwc,true,ho);
        else if (hwc==1)
            VECEXT::iirButterworth(dat,lwc,false,ho);
        else
            VECEXT::iirBpButterworth(dat,lwc,hwc,ho);
    }
    else if (type=="leaky_integration"){
        dat=std::make_shared<float1DReg> (length);
        dat->zero();
        dat->getVals()[0] = 1;
        for (int i=1; i<dat->getHyper()->getN123(); i++)
            dat->getVals()[i] = rho*dat->getVals()[i-1];
    }
    else {
        throw std::logic_error("The requested type is not implemented.\n");
    }

// adjust the phase
    if (phase=="zero"){
        dat = VECEXT::zero_phase(dat);
    }
    else if (phase=="minimum"){
        dat = VECEXT::minimum_phase(dat);
    }
    else{
        // std::clog << "The phase will be set to default.\n";
    }
        
//  set the hypercube
    axis Z;
	Z.o=(-floor(nz/2)+shift)*dz; // z origin
	Z.d=dz; // (in sec)
	Z.n=2*floor(nz/2) + 1; // number of samples in z
    std::shared_ptr<hypercube> hyper(new hypercube(Z));
    dat->setHyper(hyper);
    
// scale the wavelet
    dat->scale(amp);

// write the data to disk
    VECEXT::sepWrite(dat, output);
    
    return 0;
}