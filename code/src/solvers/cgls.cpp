#include <fstream>
#include "cgls.h"
#include "floatHyperExt.h"


using namespace SEP;

void cgls::run(const oper1D * op, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat){
    
std::clog << "Running CG solver\n";

    //data residuals
    std::shared_ptr<float1DReg> dr = dat->clone();
    dr->scale(-1);
    op->forward(true,mod,dr);

    //gradients
    std::shared_ptr<float1DReg> g (new float1DReg(mod->getHyper()));
    std::shared_ptr<float1DReg> g1 (new float1DReg(mod->getHyper()));
    op->adjoint(false,g,dr);

    //conjugate direction
    std::shared_ptr<float1DReg> p (new float1DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(dr));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary containers
    std::shared_ptr<float1DReg> d (new float1DReg(dat->getHyper()));
    std::shared_ptr<float1DReg> g_temp ;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();

    double temp;
std::clog << "Iteration = 0; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (*p->_mat)[i] = beta*(*p->_mat)[i] - (*g->_mat)[i];
        }

        //step 2: compute alpha
        op->forward(false,p,d);
        alpha = VECEXT::norm2(g) / VECEXT::norm2(d);

        //step 3: update model
        for (int i=0; i<n; i++){
            (*mod->_mat)[i] += alpha*(*p->_mat)[i];
        }

        //step 4: update data residuals
        for (int i=0; i<nd; i++){
            (*dr->_mat)[i] += alpha*(*d->_mat)[i];
        }
        
        //store the normalized data residuals
        temp = VECEXT::norm2(dr);
std::clog << "Iteration = "<<k+1<<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        op->adjoint(false,g1,dr);

        //step 6: compute beta
        beta = VECEXT::norm2(g1) / VECEXT::norm2(g);

        //copy new gradient into previous one
        g_temp = g;
        g = g1;
        g1 = g_temp;

        //step 7: iterate
        k++;
    }
    _func[0]=1;

    std::clog << "Total number of CG iterations: "<<k<<"\n";

    _res = std::dynamic_pointer_cast<float1DReg>(dr);
    _grad = std::dynamic_pointer_cast<float1DReg>(g);
}

void cgls::run(const oper2D * op, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat){
    
std::clog << "Running CG solver\n";

    //data residuals
    std::shared_ptr<float2DReg> dr = dat->clone();
    dr->scale(-1);
    op->forward(true,mod,dr);

    //gradients
    std::shared_ptr<float2DReg> g (new float2DReg(mod->getHyper()));
    std::shared_ptr<float2DReg> g1 (new float2DReg(mod->getHyper()));
    op->adjoint(false,g,dr);

    //conjugate direction
    std::shared_ptr<float2DReg> p (new float2DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(dr));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary containers
    std::shared_ptr<float2DReg> d (new float2DReg(dat->getHyper()));
    std::shared_ptr<float2DReg> g_temp ;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();

    double temp;
std::clog << "Iteration = 0; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (p->getVals())[i] = beta*(p->getVals())[i] - (g->getVals())[i];
        }

        //step 2: compute alpha
        op->forward(false,p,d);
        alpha = VECEXT::norm2(g) / VECEXT::norm2(d);

        //step 3: update model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] += alpha*(p->getVals())[i];
        }

        //step 4: update data residuals
        for (int i=0; i<nd; i++){
                (dr->getVals())[i] += alpha*(d->getVals())[i];
        }

        //store the normalized data residuals
        temp = VECEXT::norm2(dr);
std::clog << "Iteration = "<<k+1<<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        op->adjoint(false,g1,dr);

        //step 6: compute beta
        beta = VECEXT::norm2(g1) / VECEXT::norm2(g);

        //copy new gradient into previous one
        g_temp = g;
        g = g1;
        g1 = g_temp;

        //step 7: iterate
        k++;
    }

    _func[0]=1;

    std::clog << "Total number of CG iterations: "<<k<<"\n";

    _res = std::dynamic_pointer_cast<float2DReg>(dr);
    _grad = std::dynamic_pointer_cast<float2DReg>(g);
}