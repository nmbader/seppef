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
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program generate_plane_wave \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the parameters
    int nx, nz, nfx, nfz;
    float ox, oz, dx, dz, z0, vel, amp;
    bool interp;
    std::string xlabel, zlabel, interp_type;
    MB::readParam(argc, argv, "n1", nz, 1);
    MB::readParam(argc, argv, "n2", nx, 1);
    MB::readParam(argc, argv, "d1", dz, 1);
    MB::readParam(argc, argv, "d2", dx, 1);
    MB::readParam(argc, argv, "o1", oz, 0);
    MB::readParam(argc, argv, "o2", ox, 0);
    MB::readParam(argc, argv, "xlabel", xlabel, "Offset (m)");
    MB::readParam(argc, argv, "zlabel", zlabel, "Time (s)");
    MB::readParam(argc, argv, "z0", z0, 0);
    MB::readParam(argc, argv, "vel", vel, 1500);
    MB::readParam(argc, argv, "amp", amp, 1);
    MB::readParam(argc, argv, "interp", interp, false);
    MB::readParam(argc, argv, "sinc_half_lengthx", nfx, 10);
    MB::readParam(argc, argv, "sinc_half_lengthz", nfz, 10);
    MB::readParam(argc, argv, "interp_type", interp_type, "box");

// initializing the data
	axis X, Z;
	X.o=ox; // x origin (in m)
	Z.o=oz; // z origin (in m)
	X.d=dx; // x sampling (in m)
	Z.d=dz; // (in sec)
	X.n=nx; // number of samples in x
	Z.n=nz; // number of samples in z
    X.label = xlabel;
    Z.label = zlabel;
   
	std::shared_ptr<float2DReg> dat (new float2DReg(Z, X));

// compute slowness
    double s;
    if (vel >= 9999) {s = 0;}
    else {s = 1.0/vel;}
	
// filling data with a plane wave, subject to sinc interpolation if interp=true
	int iz0, min_nz, max_nz, min_nx, max_nx;
    float z;
    z0 = z0-Z.o;

    double x2, z2, ax, az, fac;
    ax = nfx*dx;
    az = nfz*dz;

    if ((interp==true) && (interp_type=="ellipse") && (nfx*nfz==0))
        throw std::logic_error("The elliptical interpolation filter length cannot be zero in either dimensions");

    if ((interp==true) && (interp_type=="box")){
        for (int ix0=0; ix0<nx; ix0++){
            z = z0 + (ix0*dx)*s;
            iz0 = z/dz;
            min_nz = std::max(0, iz0 - nfz);
            max_nz = std::min(nz, iz0 + nfz + 1);
            min_nx = std::max(0, ix0 - nfx);
            max_nx = std::min(nx, ix0 + nfx + 1);

            // interpolate with sinc(x)*sinc(z)
            for (int ix=min_nx; ix < max_nx; ix++){
                x2 = (ix0 - ix)*dx*M_PI/dx;
                for (int iz=min_nz; iz < max_nz; iz++){
                    z2 = (z - iz*dz)*M_PI/dz;
                    (*dat->_mat)[ix][iz] += amp*sinc(x2)*sinc(z2);
                }
            }
        }
    }
    else if ((interp==true) && (interp_type=="ellipse")){
        for (int ix0=0; ix0<nx; ix0++){
            z = z0 + (ix0*dx)*s;
            iz0 = z/dz;
            min_nz = std::max(0, iz0 - nfz);
            max_nz = std::min(nz, iz0 + nfz + 1);
            min_nx = std::max(0, ix0 - nfx);
            max_nx = std::min(nx, ix0 + nfx + 1);

            // interpolate with sinc(sqrt(x2+z2))
            for (int ix=min_nx; ix < max_nx; ix++){
                x2 = (ix0 - ix)*M_PI;
                x2 *=x2;
                fac = az * sqrt(1-(ix0-ix)*(ix0-ix)/(nfx*nfx));
                for (int iz=std::max(min_nz, (int)((z - fac)/dz)); iz < std::min(max_nz, (int)((z + fac)/dz)); iz++){
                    z2 = (z - iz*dz)*M_PI/dz;
                    z2 *= z2;
                    (*dat->_mat)[ix][iz] += amp*sinc(sqrt(x2 + z2));
                }
            }
        }
    }
    else{
        for (int ix=0; ix<nx; ix++){
            z = z0 + (ix*dx)*s;
            iz0 = z/dz;
            if (iz0<nz){
                (*dat->_mat)[ix][iz0] = amp;
            }
        }
    }

// write the data to disk
    VECEXT::sepWrite(dat, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program generate_plane_wave; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}