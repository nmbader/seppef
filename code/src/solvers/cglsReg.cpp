#include <fstream>
#include "cglsReg.h"
#include "floatHyperExt.h"


using namespace SEP;

void cglsReg::run(const oper1D * L, const oper1D * D, const float eps, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat, std::shared_ptr<float1DReg> mod0){
    
std::clog << "Running CG solver\n";

    //extended data residuals
    std::shared_ptr<float1DReg> drd = dat->clone();
    drd->scale(-1);
    L->forward(true,mod,drd);
    if (mod0 == nullptr){
        mod0 = std::make_shared<float1DReg>(D->getDomain());
        mod0->zero();
    }
    std::shared_ptr<float1DReg> drm (new float1DReg(D->getRange()));
    {
        std::shared_ptr<float1DReg> temp = mod->clone();
        temp->scaleAdd(mod0, 1, -1);
        D->forward(false, temp, drm);
        drm->scale(eps);
    }

    //gradients
    std::shared_ptr<float1DReg> g (new float1DReg(mod->getHyper()));
    std::shared_ptr<float1DReg> g1 (new float1DReg(mod->getHyper()));
    D->adjoint(false,g,drm);
    g->scale(eps);
    L->adjoint(true,g,drd);
    
    //conjugate direction
    std::shared_ptr<float1DReg> p (new float1DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(drd)+VECEXT::norm2(drm));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;
    
    // temporary containers
    std::shared_ptr<float1DReg> dd (new float1DReg(dat->getHyper()));
    std::shared_ptr<float1DReg> dm (new float1DReg(D->getRange()));
    std::shared_ptr<float1DReg> g_temp;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();
    unsigned int nm = drm->getHyper()->getN123();

    double temp, temp1, temp2;
    
std::clog << "Iteration = 0; d_res^2 = "<<VECEXT::norm2(drd) << "; m_res^2 = "<<VECEXT::norm2(drm)<<"; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (p->getVals())[i] = beta*(p->getVals())[i] - (g->getVals())[i];
        }
        
        //step 2: compute alpha
        L->forward(false,p,dd);
        D->forward(false,p,dm);
        dm->scale(eps);
        alpha = VECEXT::norm2(g) / (VECEXT::norm2(dd)+VECEXT::norm2(dm));
        
        //step 3: update model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] += alpha*(p->getVals())[i];
        }
        
        //step 4: update extended data residuals
        for (int i=0; i<nd; i++){
                (drd->getVals())[i] += alpha*(dd->getVals())[i];
        }
        for (int i=0; i<nm; i++){
            (drm->getVals())[i] += alpha*(dm->getVals())[i];
        }

        //store the normalized extended data residuals
        temp1 = VECEXT::norm2(drd);
        temp2 = VECEXT::norm2(drm);
        temp = temp1 + temp2;
std::clog << "Iteration = "<<k+1<<"; d_res^2 = "<<temp1 << "; m_res^2 = "<<temp2<<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        D->adjoint(false,g1,drm);
        g1->scale(eps);
        L->adjoint(true,g1,drd);

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

    D->adjoint(false,g1,drm);
    g1->scale(eps);
    _d_res = std::dynamic_pointer_cast<float1DReg>(drd);
    _m_res = std::dynamic_pointer_cast<float1DReg>(drm);
    _grad = std::dynamic_pointer_cast<float1DReg>(g);
    _m_grad = std::dynamic_pointer_cast<float1DReg>(g1);
}

void cglsReg::run(const oper2D * L, const oper2D * D, const float eps, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat, std::shared_ptr<float2DReg> mod0){
    
std::clog << "Running CG solver\n";

    //extended data residuals
    std::shared_ptr<float2DReg> drd = dat->clone();
    drd->scale(-1);
    L->forward(true,mod,drd);
    if (mod0 == nullptr){
        mod0 = std::make_shared<float2DReg>(D->getDomain());
        mod0->zero();
    }
    std::shared_ptr<float2DReg> drm (new float2DReg(D->getRange()));
    {
        std::shared_ptr<float2DReg> temp = mod->clone();
        temp->scaleAdd(mod0, 1, -1);
        D->forward(false, temp, drm);
        drm->scale(eps);
    }

    //gradients
    std::shared_ptr<float2DReg> g (new float2DReg(mod->getHyper()));
    std::shared_ptr<float2DReg> g1 (new float2DReg(mod->getHyper()));
    D->adjoint(false,g,drm);
    g->scale(eps);
    L->adjoint(true,g,drd);
    
    //conjugate direction
    std::shared_ptr<float2DReg> p (new float2DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(drd)+VECEXT::norm2(drm));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;
    
    // temporary containers
    std::shared_ptr<float2DReg> dd (new float2DReg(dat->getHyper()));
    std::shared_ptr<float2DReg> dm (new float2DReg(D->getRange()));
    std::shared_ptr<float2DReg> g_temp;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd = dat->getHyper()->getN123();
    unsigned int nm = drm->getHyper()->getN123();

    double temp, temp1, temp2;

std::clog << "Iteration = 0; d_res^2 = "<<VECEXT::norm2(drd) << "; m_res^2 = "<<VECEXT::norm2(drm)<<"; res^2 = "<<_func[0]<<"; normalized res^2 = 1\n";

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (p->getVals())[i] = beta*(p->getVals())[i] - (g->getVals())[i];
        }
        
        //step 2: compute alpha
        L->forward(false,p,dd);
        D->forward(false,p,dm);
        dm->scale(eps);
        alpha = VECEXT::norm2(g) / (VECEXT::norm2(dd)+VECEXT::norm2(dm));
        
        //step 3: update model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] += alpha*(p->getVals())[i];
        }
        
        //step 4: update extended data residuals
        for (int i=0; i<nd; i++){
                (drd->getVals())[i] += alpha*(dd->getVals())[i];
        }
        for (int i=0; i<nm; i++){
            (drm->getVals())[i] += alpha*(dm->getVals())[i];
        }

        //store the normalized extended data residuals
        temp1 = VECEXT::norm2(drd);
        temp2 = VECEXT::norm2(drm);
        temp = temp1 + temp2;
std::clog << "Iteration = "<<k+1<<"; d_res^2 = "<<temp1 << "; m_res^2 = "<<temp2<<"; res^2 = "<<temp<<"; ";
        temp = temp / _func[0];
std::clog <<"normalized res^2 = "<<temp<<"\n";
        _func.push_back(temp);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        D->adjoint(false,g1,drm);
        g1->scale(eps);
        L->adjoint(true,g1,drd);

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

    D->adjoint(false,g1,drm);
    g1->scale(eps);
    _d_res = std::dynamic_pointer_cast<float2DReg>(drd);
    _m_res = std::dynamic_pointer_cast<float2DReg>(drm);
    _grad = std::dynamic_pointer_cast<float2DReg>(g);
    _m_grad = std::dynamic_pointer_cast<float2DReg>(g1);
}