#include "oper1D.h"

void oper1D::dotProduct(){

    std::shared_ptr<float1DReg> m (new float1DReg(_domain));
    std::shared_ptr<float1DReg> mtild = m->clone();
    std::shared_ptr<float1DReg> d (new float1DReg(_range));
    std::shared_ptr<float1DReg> dtild = d->clone();

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