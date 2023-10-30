#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "floatHyperExt.h"
#include "float2DReg.h"
#include "functions.h"
#include "seplib.h"

#include "plotting.h"


using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program plot2d \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// define default ranges
    double rms = VECEXT::rms(dat);

    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);

    float xmin0, xmax0, zmin0, zmax0;
    xmin0 = X.o;
    xmax0 = X.o+X.d*(X.n-1);
    zmax0 = Z.o+Z.d*(Z.n-1);
    zmin0 = Z.o;

// read in the plot parameters
    std::string format, xlabel, zlabel, loadpath, colormap;
    float min, max, gain, xmin, xmax, zmin, zmax;
    bool zreverse;
    MB::readParam(argc, argv, "format", format, "pdf");
    MB::readParam(argc, argv, "palette_path", loadpath, "./");
    MB::readParam(argc, argv, "colormap", colormap, "BlueWhiteRed");
    MB::readParam(argc, argv, "min", min, -rms);
    MB::readParam(argc, argv, "max", max, rms);
    MB::readParam(argc, argv, "gain", gain, 1.0);
    MB::readParam(argc, argv, "xmin", xmin, xmin0);
    MB::readParam(argc, argv, "xmax", xmax, xmax0);
    MB::readParam(argc, argv, "zmin", zmin, zmin0);
    MB::readParam(argc, argv, "zmax", zmax, zmax0);
    MB::readParam(argc, argv, "xlabel", xlabel, X.label);
    MB::readParam(argc, argv, "zlabel", zlabel, Z.label);
    MB::readParam(argc, argv, "zreverse", zreverse, true);

// generate and output the data plot
    plot2D plt(dat->getHyper());
    plt.setLoadPath(loadpath);
    plt.setColormap(colormap);
    plt.setXlabel(xlabel);
    plt.setYlabel(zlabel);
    plt.setXRange(xmin, xmax);
    plt.setYRange(zmin, zmax, zreverse);
    plt.setRange(min/gain, max/gain);
    plt.setTics();
	plt.setTerminal(format);
	plt.setOutput(output);
	plt.plot(dat);

    std::clog << "Ending program plot2d\n";
    return 0;
}