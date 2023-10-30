#ifndef complexDouble6D_reg_h
#define complexDouble6D_reg_h 1
#include <complexDoubleHyper.h>
#include "boost/multi_array.hpp"


typedef boost::multi_array<std::complex<double>, 6> complexDouble6D;
namespace SEP {
/*!
A regular sampled 6-D function with complex double storage
*/
class complexDouble6DReg : public complexDoubleHyper {
 public:
  /*!
   Create a 6-D complex double vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  complexDouble6DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
 Create a 6-D complex double vector from just lengths
      \param n1,n2,n3,n4,n5,n6 Dimensions of the hypercube

 */
  complexDouble6DReg(const int n1, const int n2, const int n3, const int n4,
               const int n5, const int n6) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
    a.push_back(SEP::axis(n5));
    a.push_back(SEP::axis(n6));

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
Create a 6-D complex double vector from axes
    \param a1,a2,a3,a4,a5,a6 Axes if the hypercube

*/
  complexDouble6DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
               const SEP::axis &a4, const SEP::axis &a5, const SEP::axis &a6) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    a.push_back(a5);
    a.push_back(a6);

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 6-D complex double vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Values to fill vector with

 */
  complexDouble6DReg(std::shared_ptr<SEP::hypercube> hyper, const complexDouble6D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 6-D complex double vector from just lengths
    \param n1,n2,n3,n4,n5,n6 Dimensions of the hypercube
    \param vals Values to fill vector with

*/
  complexDouble6DReg(const int n1, const int n2, const int n3, const int n4,
               const int n5, const int n6, complexDouble6D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
    a.push_back(SEP::axis(n5));
    a.push_back(SEP::axis(n6));

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
  Create a 6-D complex double vector from axes
       \param a1,a2,a3,a4,a5,a6 Axes if the hypercube
       \param vals Values to fill vector with
  */
  complexDouble6DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
               const SEP::axis &a4, const SEP::axis &a5, const SEP::axis &a6,
               complexDouble6D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    a.push_back(a5);
    a.push_back(a6);

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
    Allocate data for vector
    */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complexDouble6D(
        boost::extents[ns[5]][ns[4]][ns[3]][ns[2]][ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
Make a copy of the vector
*/
  std::shared_ptr<complexDouble6DReg> clone() const;
  /*!
   Make a copy of the vector space
   */
  std::shared_ptr<complexDouble6DReg> cloneSpace() const;
  /*!
   Deallocate storage for vector, turn into vector space
    */
  virtual void cleanMemory() {
    _mat = 0;
    setSpace();
  }
  /*!
                                Return a subset of the vector
                                \param nw,fw,jw Windowing parameters
                                */
  std::shared_ptr<complexDouble6DReg> window(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw) const;
  std::shared_ptr<complexDouble6D> _mat;  ///< Storage for vector

 private:
  /*!
  Initialize vector without providing data (zero)
  \param hyp Hypercube describing dataset
*/
  void initNoData(std::shared_ptr<SEP::hypercube> hyper);
  /*!
Initialize with data
\param hyper Hypercube describing space
\param vals Data to copy in
*/
  void initData(std::shared_ptr<SEP::hypercube> hyper, const complexDouble6D &vals);
};
}  // namespace SEP
#endif
