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
#include "float2DReg.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;

double xlogx(double x){
    assert (x >= 0);
    if (x==0) return 0;
    else return x*log2(x);
}


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program stat \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// list of available statistics
std::vector<std::string> list = {"xstack", "zstack", "xmax", "zmax", "xmin", "zmin", "xmaxabs", "zmaxabs", "xrms",
                                "zrms","pdf"};

// read in the statistics type
    std::string type;
    bool normalize, logabs;
    MB::readParam(argc, argv, "type", type, "zrms");
    MB::readParam(argc, argv, "normalize", normalize, false);

// compute the requested statistics
    axis Z = dat->getHyper()->getAxis(1);
    axis X = dat->getHyper()->getAxis(2);

    if (type == "zrms"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(X));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[ix] += (*dat->_mat)[ix][iz]*(*dat->_mat)[ix][iz];
            }
            val[ix] = normalize*sqrt(val[ix]/Z.n) + (1-normalize)*sqrt(val[ix]);
        }

        // write the statistics vector to disk
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "xrms"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(Z));
        out0->zero();
        float* val = out0->getVals();
        for (int iz=0; iz<Z.n; iz++){
            for (int ix=0; ix<X.n; ix++){
                val[iz] += (*dat->_mat)[ix][iz]*(*dat->_mat)[ix][iz];
            }
            val[iz] = normalize*sqrt(val[iz]/X.n) + (1-normalize)*sqrt(val[iz]);
        }

        VECEXT::sepWrite(out0, output);
    }

    else if (type == "zmax"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(X));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[ix] = std::max(val[ix], (*dat->_mat)[ix][iz]);
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "zmin"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(X));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[ix] = std::min(val[ix], (*dat->_mat)[ix][iz]);
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "zmaxabs"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(X));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[ix] = std::max(val[ix], std::abs((*dat->_mat)[ix][iz]));
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "xmax"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(Z));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[iz] = std::max(val[iz], (*dat->_mat)[ix][iz]);
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "xmin"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(Z));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[iz] = std::min(val[iz], (*dat->_mat)[ix][iz]);
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "xmaxabs"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(Z));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[iz] = std::max(val[iz], std::abs((*dat->_mat)[ix][iz]));
            }
        }
        VECEXT::sepWrite(out0, output);
    }

    else if (type == "xstack"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(Z));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[iz] += (*dat->_mat)[ix][iz];
            }
        }
        if (normalize==true){
            out0->scale(1.0/X.n);
        }
        VECEXT::sepWrite(out0, output);
        
    }

    else if (type == "zstack"){
        std::shared_ptr<float1DReg> out0 (new float1DReg(X));
        out0->zero();
        float* val = out0->getVals();
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                val[ix] += (*dat->_mat)[ix][iz];
            }
        }
        if (normalize==true){
            out0->scale(1.0/Z.n);
        }
        VECEXT::sepWrite(out0, output);
        
    }

    else if (type == "pdf"){
        assert (dat->min() >= 0);
        dat->scale(1.0/VECEXT::sum(dat)); // transform positive function to a PDF

        double xmean, zmean, xstddev, zstddev, whiteness, cte, stddev, flatness, sum;
        xmean=0; zmean=0; xstddev=0; zstddev=0; whiteness=0; stddev = 0; flatness = 0; sum = 0;

        cte = 1.0 / (X.n * Z.n); // equivalent uniform PDF

        // compute the expected value and the whiteness
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                xmean += (*dat->_mat)[ix][iz] * (ix*X.d - X.o);
                zmean += (*dat->_mat)[ix][iz] * (iz*Z.d - Z.o);
                whiteness += std::min(cte, (double)(*dat->_mat)[ix][iz]);
                stddev += ((double)(*dat->_mat)[ix][iz] - cte)*((double)(*dat->_mat)[ix][iz] - cte);
                flatness += xlogx((*dat->_mat)[ix][iz]);
                sum += (*dat->_mat)[ix][iz];
            }
        }
        stddev = sqrt(stddev/(X.n * Z.n));
        flatness /= log2(X.n * Z.n);
        flatness = pow(2,-flatness) - 1;

        // compute the std deviations
        for (int ix=0; ix<X.n; ix++){
            for (int iz=0; iz<Z.n; iz++){
                xstddev += (*dat->_mat)[ix][iz] * pow((ix*X.d - X.o - xmean),2);
                zstddev += (*dat->_mat)[ix][iz] * pow((iz*Z.d - Z.o - zmean),2);
            }
        }

std::clog << "PDF X mean = " << xmean << "\n";
std::clog << "PDF Z mean = " << zmean << "\n";
std::clog << "PDF X stddev = " << xstddev << "\n";
std::clog << "PDF Z stddev = " << zstddev << "\n";
std::clog << "PDF whiteness = " << whiteness << "\n";
std::clog << "PDF stddev = " << stddev << "\n";
std::clog << "PDF flatness = " << flatness << "\n";
std::clog << "PDF sum = " << sum << "\n";

    }

    else{
        throw std::logic_error("The requested statistics is not implemented.\n");
    }

	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program stat; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}