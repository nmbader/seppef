#include "oper2D.h"
#include "floatHyperExt.h"

void oper2D::dotProduct(){

    std::shared_ptr<float2DReg> m (new float2DReg(_domain));
    std::shared_ptr<float2DReg> mtild = m->clone();
    std::shared_ptr<float2DReg> d (new float2DReg(_range));
    std::shared_ptr<float2DReg> dtild = d->clone();

    m->random();
    d->random();

    this->forward(false, m, dtild); 
    this->adjoint(false, mtild, d);
    
    double sum1 = d->dot(dtild);
    double sum2 = m->dot(mtild);

    std::clog << "Dot product with add=false:\n";
    std::clog << "sum1 = " << sum1 << std::endl;
    std::clog << "sum2 = " << sum2 << std::endl;
    std::clog << "diff (x1000) = " << 1000*(sum1 - sum2) << std::endl;

    this->forward(true, m, dtild); 
    this->adjoint(true, mtild, d);

    sum1 = d->dot(dtild);
    sum2 = m->dot(mtild);

    std::clog << "Dot product with add=true:\n";
    std::clog << "sum1 = " << sum1 << std::endl;
    std::clog << "sum2 = " << sum2 << std::endl;
    std::clog << "diff (x1000) = " << 1000*(sum1 - sum2) << std::endl;
}