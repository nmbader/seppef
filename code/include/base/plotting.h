#ifndef plotting_H
#define plotting_H

//#include <boost/tuple/tuple.hpp>
//#include <boost/foreach.hpp>

#include <chrono>
#include <thread>
#include <stdexcept>
#include "gnuplot-iostream.h"
#include "float1DReg.h"
#include "float2DReg.h"
#include "hypercube.h"

using namespace SEP;

class plot2D {

private:
    Gnuplot _gp;
    std::shared_ptr<hypercube> _sh;
    std::string _terminal, _loadPath;

public:
    // constructor
    plot2D(){}

    // destructor
    ~plot2D(){}

    // constructor from float2DReg
    plot2D(const std::shared_ptr<float2DReg> dat){
        _sh = dat->getHyper()->clone();
        _terminal = "x11";
        _loadPath = "../misc/gnuplot_palettes/";
        this->setDefault();
    }

    // constructor from hypercube
    plot2D(const std::shared_ptr<hypercube> sh){
        _sh = sh->clone();
        this->setDefault();
    }

    // set default parameters
    void setDefault(){
        _gp << "set xrange ["<<_sh->getAxis(2).o<<":"<<_sh->getAxis(2).o+_sh->getAxis(2).d*(_sh->getAxis(2).n-1)<<"]\n";
        _gp << "set yrange ["<<_sh->getAxis(1).o+_sh->getAxis(1).d*(_sh->getAxis(1).n-1)<<":"<<_sh->getAxis(1).o<<"]\n";
        setXlabel("Offset (m)");
        setYlabel("Time (sec)");
    }

    // set load path for the palettes
    void setLoadPath(std::string loadPath){
        _loadPath = loadPath;
    }

    // set title
    void setTitle(std::string title){
        _gp << "set title \""<< title << "\" \n";
    }

    // set colormap
    void setColormap(std::string colormap){
        if (colormap == "default"){
            _gp << "set palette defined\n";
        }
        else {
        _gp << "set loadpath \""<< _loadPath <<"\"\n";
        _gp << "load '"<<colormap<<".pal'\n";
        }
    }

    // set xrange
    void setXRange(float min, float max, const bool reverse = false){
        std::string rev = "noreverse";
        if (reverse == true) rev = "reverse";
        _gp << "set xrange [" << min << ":" << max << "] " << rev << "\n";
    }

    // set yrange
    void setYRange(float min, float max, const bool reverse = true){
        std::string rev = "reverse";
        if (reverse == false) rev = "noreverse";
        _gp << "set yrange [" << min << ":" << max << "] " << rev << "\n";
    }

    // set colormap range
    void setRange(float min, float max){
        _gp << "set cbrange [" << min << ":" << max << "]\n";
    }

    // set X axis label
    void setXlabel(std::string label){
        _gp << "set xlabel \""<< label <<"\" \n";
    }

    // set Y axis label
    void setYlabel(std::string label){
        _gp << "set ylabel \""<< label <<"\" \n";
    }

    // set tics
    void setTics(){
        _gp << "unset tics \n";
        _gp << "set xtics border nomirror out \n";
        _gp << "set ytics border nomirror out \n";
        _gp << "set cbtics border nomirror out \n";
    }

    // set the terminal for gnuplot
    void setTerminal(std::string term = "x11"){
        _gp << "set terminal "<< term <<"\n";
        _terminal = term;
    }

    // set the output name and directory
    void setOutput(std::string name = "plot2D"){
        _gp << "set output '"<< name <<"'\n";
    }

    // plot the image of a float2DReg
    void plot(const std::shared_ptr<float2DReg> vec){
        _gp << "plot '-' matrix using ($2*"<<_sh->getAxis(2).d<<"+"<<_sh->getAxis(2).o<<"):($1*"<<_sh->getAxis(1).d<<"+"<<_sh->getAxis(1).o<<"):3 with image notitle\n";
        unsigned int nx = vec->getHyper()->getAxis(2).n;
        unsigned int nz = vec->getHyper()->getAxis(1).n;
        std::vector<std::vector<float> > v(nx, std::vector<float>(nz));
        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                v[ix][iz] = (*vec->_mat)[ix][iz];
            }
        }
        _gp.send1d(v);
        //this_thread::sleep_for(chrono::milliseconds(200));
    }

};


class plot1D {

private:
    Gnuplot _gp;
    std::shared_ptr<hypercube> _sh;
    std::string _terminal;

public:
    // constructor
    plot1D(){}

    // destructor
    ~plot1D(){}

    // constructor from float1DReg
    plot1D(const std::shared_ptr<float1DReg> dat){
        _sh = dat->getHyper()->clone();
        _terminal = "x11";
        this->setDefault();
    }

    // constructor from hypercube
    plot1D(const std::shared_ptr<hypercube> sh){
        _sh = sh->clone();
        this->setDefault();
    }

    // set default parameters
    void setDefault(){
        _gp << "set xrange ["<<_sh->getAxis(1).o<<":"<<_sh->getAxis(1).o+_sh->getAxis(1).d*(_sh->getAxis(1).n-1)<<"]\n";
        setXlabel("Offset (m)");
        setXlabel("Value");
    }

    // set title
    void setTitle(std::string title){
        _gp << "set title \""<< title << "\" \n";
    }

    // set xrange
    void setXRange(float min, float max, const bool reverse = false){
        std::string rev = "noreverse";
        if (reverse == true) rev = "reverse";
        _gp << "set xrange [" << min << ":" << max << "] " << rev << "\n";
    }

    // set yrange
    void setYRange(float min, float max, const bool reverse = false){
        std::string rev = "noreverse";
        if (reverse == true) rev = "reverse";
        _gp << "set yrange [" << min << ":" << max << "] " << rev << "\n";
    }

    // set X axis label
    void setXlabel(std::string label){
        _gp << "set xlabel \""<< label <<"\" \n";
    }

    // set Y axis label
    void setYlabel(std::string label){
        _gp << "set ylabel \""<< label <<"\" \n";
    }

    // set X logscale
    void setXlogscale(){
        _gp << "set logscale x\n";
    }

    // set Y logscale
    void setYlogscale(){
        _gp << "set logscale y\n";
    }

    // set tics
    void setTics(){
        _gp << "unset tics \n";
        _gp << "set xtics border nomirror out \n";
        _gp << "set ytics border nomirror out \n";
    }

    // set the terminal for gnuplot
    void setTerminal(std::string term = "x11"){
        _gp << "set terminal "<< term <<"\n";
        _terminal = term;
    }

    // set the output name and directory
    void setOutput(std::string name = "plot1D"){
        _gp << "set output '"<< name <<"'\n";
    }

    // plot a float1DReg
    void plot(const std::shared_ptr<float1DReg> vec){
        _gp << "plot '-' using ($0*"<<_sh->getAxis(1).d<<"+"<<_sh->getAxis(1).o<<"):1 with linespoints notitle\n";
        unsigned int nz = vec->getHyper()->getAxis(1).n;
        std::vector<float> v(nz);
        for (int iz=0; iz<nz; iz++){
            v[iz] = (*vec->_mat)[iz];
        }
        _gp.send1d(v);
    }

    // plot an x vector from a float2DReg
    void plotx(const std::shared_ptr<float2DReg> vec, unsigned int ix = 0){
        if (ix >= vec->getHyper()->getAxis(2).n)
            throw std::logic_error("Index ix out of range");

        _gp << "plot '-' using ($0*"<<_sh->getAxis(1).d<<"+"<<_sh->getAxis(1).o<<"):1 with linespoints notitle\n";
        unsigned int nz = vec->getHyper()->getAxis(1).n;
        std::vector<float> v(nz);
        for (int iz=0; iz<nz; iz++){
            v[iz] = (*vec->_mat)[ix][iz];
        }
        _gp.send1d(v);
    }

    // plot a z vector from a float2DReg
    void plotz(const std::shared_ptr<float2DReg> vec, unsigned int iz = 0){
        if (iz >= vec->getHyper()->getAxis(1).n)
            throw std::logic_error("Index iz out of range");

        _gp << "set xrange ["<<_sh->getAxis(2).d<<":"<<_sh->getAxis(2).o+_sh->getAxis(2).d*(_sh->getAxis(2).n-1)<<"]\n";
        _gp << "plot '-' using ($0*"<<_sh->getAxis(2).d<<"+"<<_sh->getAxis(2).o<<"):1 with linespoints notitle\n";

        unsigned int nx = vec->getHyper()->getAxis(2).n;
        std::vector<float> v(nx);
        for (int ix=0; ix<nx; ix++){
            v[ix] = (*vec->_mat)[ix][iz];
        }
        _gp.send1d(v);
    }

};

#endif