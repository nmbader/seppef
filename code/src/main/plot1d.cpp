#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "floatHyperExt.h"
#include "float1DReg.h"
#include "functions.h"
#include "seplib.h"

#include "plotting.h"


using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program plot1d \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float1DReg> dat = VECEXT::sepRead1D();

// define default ranges
    double min = dat->min();
    double max = dat->max();

    axis X = dat->getHyper()->getAxis(1);

    float xmin0, xmax0;
    xmin0 = X.o;
    xmax0 = X.o+X.d*(X.n-1);

// read in the plot parameters
    std::string format, xlabel, zlabel;
    float xmin, xmax, zmin, zmax;
    bool zreverse, xlogscale, zlogscale;
    MB::readParam(argc, argv, "format", format, "pdf");
    MB::readParam(argc, argv, "xmin", xmin, xmin0);
    MB::readParam(argc, argv, "xmax", xmax, xmax0);
    MB::readParam(argc, argv, "zmin", zmin, min);
    MB::readParam(argc, argv, "zmax", zmax, max);
    MB::readParam(argc, argv, "xlabel", xlabel, X.label);
    MB::readParam(argc, argv, "zlabel", zlabel, "Value");
    MB::readParam(argc, argv, "zreverse", zreverse, false);
    MB::readParam(argc, argv, "xlogscale", xlogscale, false);
    MB::readParam(argc, argv, "zlogscale", zlogscale, false);

// generate and output the data plot
    plot1D plt(dat->getHyper());
    plt.setXlabel(xlabel);
    plt.setYlabel(zlabel);
    plt.setXRange(xmin, xmax);
    plt.setYRange(zmin, zmax, zreverse);
    if (xlogscale==true)
        plt.setXlogscale();
    if (zlogscale==true)
        plt.setYlogscale();
    plt.setTics();
	plt.setTerminal(format);
	plt.setOutput(output);
	plt.plot(dat);

    std::clog << "Ending program plot1d\n";
    return 0;
}