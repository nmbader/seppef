#ifndef complexDouble3D_reg_h
#define complexDouble3D_reg_h 1
#include <complexDoubleHyper.h>
#include "boost/multi_array.hpp"
typedef boost::multi_array<std::complex<double>, 3> complexDouble3D;
//! Namespace for SEP code
namespace SEP {
/*!
A regular sampled 3-D function with complex double storage
*/
class complexDouble3DReg : public complexDoubleHyper {
 public:
  /*!
   Create a 3-D complex double vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  complexDouble3DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 3-D complex double vector from just lengths
     \param n1,n2,n3 Dimensions of the hypercube

*/
  complexDouble3DReg(const int n1, const int n2, const int n3) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 3-D complex double vector from axes
      \param a1,a2,a3 Axes if the hypercube

 */
  complexDouble3DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 3-D complex double vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Vaule to file vector with
 */
  complexDouble3DReg(std::shared_ptr<SEP::hypercube> hyper, const complexDouble3D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 3-D complex double vector from just lengths
  \param n1,n2,n3 Dimensions of the hypercube
  \param vals Vaule to file vector with

*/
  complexDouble3DReg(int n1, int n2, int n3, complexDouble3D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
Create a 3-D complex double vector from axes
    \param a1,a2,a3 Axes if the hypercube
    \param vals Vaule to file vector with

*/
  complexDouble3DReg(SEP::axis &a1, SEP::axis &a2, SEP::axis &a3, complexDouble3D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
  Allocate data for vector
  */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complexDouble3D(boost::extents[ns[2]][ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
   Make a copy of the vector
   */
  std::shared_ptr<complexDouble3DReg> clone() const;
  /*!
Make a copy of the vector space
*/
  std::shared_ptr<complexDouble3DReg> cloneSpace() const;
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
  std::shared_ptr<complexDouble3DReg> window(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw) const;
  std::shared_ptr<complexDouble3D> _mat;  ///< Storage for vector

 protected:
  /*!
 Initialize without data
 \param hyper Hypercube describing space
*/
  void initNoData(std::shared_ptr<SEP::hypercube> hyper);
  /*!
  Initialize with data
  \param hyper Hypercube describing space
  \param vals Data to copy in
*/

  void initData(std::shared_ptr<SEP::hypercube> hyper, const complexDouble3D &vals);
  /*!
Initialize with data
\param hyper Hypercube describing space
\param vals Data to copy in
*/

  void initData(std::shared_ptr<SEP::hypercube> hyper,
                std::shared_ptr<complexDouble3D> vals);
};
}  // namespace SEP
#endif
