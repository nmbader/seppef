#include <fstream>
#include "sdls.h"
#include "floatHyperExt.h"

using namespace SEP;

void sdls::run(const oper1D * op, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat){
    
std::clog << "Running SD solver\n";

    //data residuals
    std::shared_ptr<float1DReg> res = dat->clone();
    res->scale(-1);
    op->forward(true,mod,res);

    //gradient
    std::shared_ptr<float1DReg> grad (new float1DReg(mod->getHyper()));
    op->adjoint(false,grad,res);

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(res));

    //step length
    double alpha = 0;

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary container
    std::shared_ptr<float1DReg> d (new float1DReg(dat->getHyper()));

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();

    double temp;

std::clog << "Iteration = 0; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the SD loop consisting of 4 steps
    while (k<_niter & rate>_threshold){

        //step 1: alpha = grad^2 / (Op*grad)^2
        op->forward(false,grad,d);
        alpha = VECEXT::norm2(grad) / VECEXT::norm2(d);

        //step 2: update the model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] -= alpha*(grad->getVals())[i];
        }

        //step 3: update data residuals
        op->forward(false,mod,res);
        for (int i=0; i<nd; i++){
            (res->getVals())[i] -= (dat->getVals())[i];
        }

        //step4: compute the gradient
        op->adjoint(false,grad,res);

        //store the normalized data residuals
        temp = VECEXT::norm2(res);
std::clog << "Iteration = "<<k+1<<"; Step Length = "<< alpha <<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //iterate
        k++;
    }
    _func[0]=1;

    std::clog << "Total number of SD iterations: "<<k<<"\n";

    _res = std::dynamic_pointer_cast<float1DReg>(res);
    _grad = std::dynamic_pointer_cast<float1DReg>(grad);
}

void sdls::run(const oper2D * op, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat){
    
std::clog << "Running SD solver\n";

    //data residuals
    std::shared_ptr<float2DReg> res = dat->clone();
    res->scale(-1);
    op->forward(true,mod,res);

    //gradient
    std::shared_ptr<float2DReg> grad (new float2DReg(mod->getHyper()));
    op->adjoint(false,grad,res);

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(res));

    //step length
    double alpha = 0;

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary container
    std::shared_ptr<float2DReg> d (new float2DReg(dat->getHyper()));

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();

    double temp;

std::clog << "Iteration = 0; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the SD loop consisting of 4 steps
    while (k<_niter & rate>_threshold){

        //step 1: alpha = grad^2 / (Op*grad)^2
        op->forward(false,grad,d);
        alpha = VECEXT::norm2(grad) / VECEXT::norm2(d);

        //step 2: update the model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] -= alpha*(grad->getVals())[i];
        }

        //step 3: update data residuals
        op->forward(false,mod,res);
        for (int i=0; i<nd; i++){
            (res->getVals())[i] -= (dat->getVals())[i];
        }

        //step4: compute the gradient
        op->adjoint(false,grad,res);

        //store the normalized data residuals
        temp = VECEXT::norm2(res);
std::clog << "Iteration = "<<k+1<<"; Step Length = "<< alpha <<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //iterate
        k++;
    }
    _func[0]=1;

    std::clog << "Total number of SD iterations: "<<k<<"\n";

    _res = std::dynamic_pointer_cast<float2DReg>(res);
    _grad = std::dynamic_pointer_cast<float2DReg>(grad);
}