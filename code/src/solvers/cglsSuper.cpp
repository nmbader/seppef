#include "cglsSuper.h"
#include "matrixMult.h"
#include "identity1D.h"
#include "floatHyperExt.h"

using namespace SEP;

void cglsSuper::run(const oper1D * L1, const oper1D * L2, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat1, const std::shared_ptr<float1DReg> dat2){
    //initialize variables

    //extended data residuals
    std::shared_ptr<float1DReg> drd1 = dat1->clone();
    drd1->scale(-1);
    L1->forward(true,mod,drd1);
    std::shared_ptr<float1DReg> drd2 = dat2->clone();
    drd2->scale(-1);
    L2->forward(true,mod,drd2);

    //gradients
    std::shared_ptr<float1DReg> g (new float1DReg(mod->getHyper()));
    std::shared_ptr<float1DReg> g1 (new float1DReg(mod->getHyper()));
    L1->adjoint(false,g,drd1);
    L2->adjoint(true,g,drd2);
    
    //conjugate direction
    std::shared_ptr<float1DReg> p (new float1DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(drd1)+VECEXT::norm2(drd2));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary containers
    std::shared_ptr<float1DReg> dd1 (new float1DReg(dat1->getHyper()));
    std::shared_ptr<float1DReg> dd2 (new float1DReg(dat2->getHyper()));
    std::shared_ptr<float1DReg> g_temp;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd1 = dat1->getHyper()->getN123();
    unsigned int nd2 = dat2->getHyper()->getN123();

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (p->getVals())[i] = beta*(p->getVals())[i] - (g->getVals())[i];
        }
        
        //step 2: compute alpha
        L1->forward(false,p,dd1);
        L2->forward(false,p,dd2);

        alpha = VECEXT::norm2(g) / (VECEXT::norm2(dd1)+VECEXT::norm2(dd2));
        
        //step 3: update model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] += alpha*(p->getVals())[i];
        }
        
        //step 4: update extended data residuals
        for (int i=0; i<nd1; i++){
                (drd1->getVals())[i] += alpha*(dd1->getVals())[i];
        }
        for (int i=0; i<nd2; i++){
                (drd2->getVals())[i] += alpha*(dd2->getVals())[i];
        }

        //store the normalized extended data residuals
        _func.push_back((VECEXT::norm2(drd1)+VECEXT::norm2(drd2))/_func[0]);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        L1->adjoint(false,g1,drd1);
        L2->adjoint(true,g1,drd2);

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
}

void cglsSuper::run(const oper2D * L1, const oper2D * L2, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat1, const std::shared_ptr<float2DReg> dat2){
    //initialize variables

    //extended data residuals
    std::shared_ptr<float2DReg> drd1 = dat1->clone();
    drd1->scale(-1);
    L1->forward(true,mod,drd1);
    std::shared_ptr<float2DReg> drd2 = dat2->clone();
    drd2->scale(-1);
    L2->forward(true,mod,drd2);

    //gradients
    std::shared_ptr<float2DReg> g (new float2DReg(mod->getHyper()));
    std::shared_ptr<float2DReg> g1 (new float2DReg(mod->getHyper()));
    L1->adjoint(false,g,drd1);
    L2->adjoint(true,g,drd2);
    
    //conjugate direction
    std::shared_ptr<float2DReg> p (new float2DReg(mod->getHyper()));

    //directions weights
    float alpha = 0;
    float beta = 0;

    //objective function
    _func.clear();
    _func.push_back(VECEXT::norm2(drd1)+VECEXT::norm2(drd2));

    //iteration number
    unsigned int k = 0;

    //convergence rate
    float rate = 1;

    // temporary containers
    std::shared_ptr<float2DReg> dd1 (new float2DReg(dat1->getHyper()));
    std::shared_ptr<float2DReg> dd2 (new float2DReg(dat2->getHyper()));
    std::shared_ptr<float2DReg> g_temp;

    unsigned int n = mod->getHyper()->getN123();
    unsigned int nd1 = dat1->getHyper()->getN123();
    unsigned int nd2 = dat2->getHyper()->getN123();

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        for (int i=0; i<n; i++){
            (p->getVals())[i] = beta*(p->getVals())[i] - (g->getVals())[i];
        }
        
        //step 2: compute alpha
        L1->forward(false,p,dd1);
        L2->forward(false,p,dd2);

        alpha = VECEXT::norm2(g) / (VECEXT::norm2(dd1)+VECEXT::norm2(dd2));
        
        //step 3: update model
        for (int i=0; i<n; i++){
            (mod->getVals())[i] += alpha*(p->getVals())[i];
        }
        
        //step 4: update extended data residuals
        for (int i=0; i<nd1; i++){
                (drd1->getVals())[i] += alpha*(dd1->getVals())[i];
        }
        for (int i=0; i<nd2; i++){
                (drd2->getVals())[i] += alpha*(dd2->getVals())[i];
        }

        //store the normalized extended data residuals
        _func.push_back((VECEXT::norm2(drd1)+VECEXT::norm2(drd2))/_func[0]);

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
            if((_func[k]<1e-07) || (_func[k]-_func[k+1]<1e-07)) rate = 0;
        }

        //step 5: update gradient
        L1->adjoint(false,g1,drd1);
        L2->adjoint(true,g1,drd2);

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
}

void cglsSuper::test(){

    matrixMult mult1(50,200);
    matrixMult mult2(50,150);
    std::shared_ptr<float1DReg> mod (new float1DReg(50));
    std::shared_ptr<float1DReg> mod0 (new float1DReg(50));
    std::shared_ptr<float1DReg> dat1 (new float1DReg(200));
    std::shared_ptr<float1DReg> dat2 (new float1DReg(150));

    mult1.random(-1,1);
    mult2.random(-1,1);
    VECEXT::random(mod,-10,10);
    mult1.forward(false,mod,dat1);
    mult2.forward(false,mod,dat2);

    _niter = 50;
    _threshold = 0.0;

    this->run(&mult1, &mult2, mod0, dat1, dat2);

    // L2 squared error of estimated model
    double trueNorm, estimatedNorm, dataNorm, error;
    dataNorm = VECEXT::norm2(dat1) + VECEXT::norm2(dat2);
    trueNorm = VECEXT::norm2(mod);
    estimatedNorm = VECEXT::norm2(mod0);
    mod->scale(-1);
    mod0->add(mod);
    error = VECEXT::norm2(mod0);

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "estimated model norm squared = " << estimatedNorm << std::endl;
    std::clog << "model error squared = " << error << std::endl;
}